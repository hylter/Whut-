#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#define CLIGHT      299792458.0 
// GPS Frequencies
#define FREQ_L1     1.57542e9    
#define FREQ_L2     1.22760e9    
// Galileo Frequencies
#define FREQ_E1     1.57542e9  
#define FREQ_E5a    1.17645e9  
// [重点修正] BDS Frequencies
// 依据您的指示：B1 优先指代 1575.42 MHz (B1C)
#define FREQ_B1I    1.561098e9 // B1I (D2I)
#define FREQ_B1C    1.57542e9  // B1C (D1X/D1I) - 1575.42 MHz
#define FREQ_B3I    1.268520e9 // B3I (D6I)     - 1268.52 MHz

#define GM_GPS      3.986005e14     
#define GM_GAL      3.986004418e14  
#define GM_BDS      3.986004418e14  
#define OMEGA_E_BDS 7.292115e-5     

#define OMEGA_E     7.2921151467e-5  
#define PI          3.1415926535897932
#define MAX_SAT     64       
#define MAX_SYS     3        
#define MAX_EPOCH   2880     
#define MAX_NAV_PER_SAT 200   

#define SYS_GPS 0   
#define SYS_GAL 1   
#define SYS_BDS 2   

typedef struct {
    time_t time;
    double sec;
} gtime_t;

typedef struct {
    int prn;
    int sys;
    gtime_t toc;
    gtime_t toe;
    double toe_sec;
    double sqrtA, e, i0, OMEGA0, omega, M0, dn, OMEGAdot, IDOT;
    double cuc, cus, crc, crs, cic, cis, a0, a1, a2;
} nav_t;

typedef struct {
    int prn;
    int sys;
    double D1, D2;
} sat_obs_t;

typedef struct {
    gtime_t time;
    int sat_num;
    sat_obs_t sats[MAX_SAT * MAX_SYS];
} epoch_obs_t;

typedef struct {
    double x, y, z;
} vec3_t;

//  全局变量
nav_t nav_data[MAX_SAT * MAX_SYS + 1][MAX_NAV_PER_SAT];
int nav_count[MAX_SAT * MAX_SYS + 1] = { 0 };
epoch_obs_t obs_data[MAX_EPOCH];
int obs_count = 0;
vec3_t STA_POS_XYZ = { 0 };

int idx_G_D1 = -1, idx_G_D2 = -1;
int idx_E_D1 = -1, idx_E_D2 = -1;
int idx_C_D1 = -1, idx_C_D2 = -1;

// 动态记录 BDS B1 频率
double cur_bds_f1 = FREQ_B1C; // 默认为 B1C (1575.42)

void str_replace_d_e(char* str) {
    for (; *str; str++) if (*str == 'D' || *str == 'd') *str = 'E';
}

gtime_t epoch2time(const double* ep) {
    const int doy[] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
    gtime_t time = { 0 };
    int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];
    if (year < 1970) year += 2000;
    days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
    sec = (int)ep[3] * 3600 + (int)ep[4] * 60 + (int)ep[5];
    time.time = (time_t)days * 86400 + sec;
    time.sec = ep[5] - (int)ep[5];
    return time;
}

double time2gpst(gtime_t t, int* week) {
    double ep_0[6] = { 1980, 1, 6, 0, 0, 0 };
    gtime_t t0 = epoch2time(ep_0);
    time_t sec = t.time - t0.time;
    int w = (int)(sec / (86400 * 7));
    if (week) *week = w;
    return (double)(sec - w * 86400 * 7) + t.sec;
}

double timediff(gtime_t t1, gtime_t t2) {
    return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}

int matinv4(double* A, double* invA) {
    int i, j, k;
    double B[4][8], ratio;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) B[i][j] = A[i * 4 + j];
        for (j = 4; j < 8; j++) B[i][j] = (i == j - 4) ? 1.0 : 0.0;
    }
    for (i = 0; i < 4; i++) {
        if (fabs(B[i][i]) < 1e-12) return 0;
        for (j = 0; j < 4; j++) {
            if (i != j) {
                ratio = B[j][i] / B[i][i];
                for (k = 0; k < 8; k++) B[j][k] -= ratio * B[i][k];
            }
        }
    }
    for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) invA[i * 4 + j] = B[i][j + 4] / B[i][i];
    return 1;
}

void xyz2enu(vec3_t ref, vec3_t diff, double* E, double* N, double* U) {
    double r = sqrt(ref.x * ref.x + ref.y * ref.y);
    double lat = atan2(ref.z, r);
    double lon = atan2(ref.y, ref.x);
    double sin_p = sin(lat), cos_p = cos(lat), sin_l = sin(lon), cos_l = cos(lon);
    *E = -sin_l * diff.x + cos_l * diff.y;
    *N = -sin_p * cos_l * diff.x - sin_p * sin_l * diff.y + cos_p * diff.z;
    *U = cos_p * cos_l * diff.x + cos_p * sin_l * diff.y + sin_p * diff.z;
}

int get_idx(int prn, int sys) {
    if (sys == SYS_GAL) return prn + MAX_SAT;
    if (sys == SYS_BDS) return prn + MAX_SAT * 2;
    return prn;
}

// 文件读取模块
void read_nav(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) { printf("错误: 无法打开 NAV 文件: %s\n", filename); return; }
    char line[1024];
    int prn = 0, sys = SYS_GPS;
    int nav_total = 0;

    while (fgets(line, sizeof(line), fp)) if (strstr(line, "END OF HEADER")) break;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == 'G') { sys = SYS_GPS; prn = atoi(line + 1); }
        else if (line[0] == 'E') { sys = SYS_GAL; prn = atoi(line + 1); }
        else if (line[0] == 'C') { sys = SYS_BDS; prn = atoi(line + 1); }
        else if (isdigit(line[0])) { sys = SYS_GPS; prn = atoi(line); }
        else continue;

        if (prn <= 0 || prn > MAX_SAT) continue;
        int idx = get_idx(prn, sys);

        double ep[6];
        str_replace_d_e(line);
        char* p_time = (line[0] == 'G' || line[0] == 'E' || line[0] == 'C') ? line + 3 : line;
        if (sscanf(p_time, "%lf %lf %lf %lf %lf %lf", &ep[0], &ep[1], &ep[2], &ep[3], &ep[4], &ep[5]) < 6) continue;

        int count = nav_count[idx];
        if (count >= MAX_NAV_PER_SAT) {
            for (int k = 0; k < 7; k++) if (!fgets(line, sizeof(line), fp)) break;
            continue;
        }
        nav_t* nav = &nav_data[idx][count];
        nav->prn = prn;
        nav->sys = sys;
        nav->toc = epoch2time(ep);
        nav->toe = nav->toc;

        char* p = line + 23;
        nav->a0 = strtod(p, &p); nav->a1 = strtod(p, &p); nav->a2 = strtod(p, &p);

        if (!fgets(line, sizeof(line), fp)) break; str_replace_d_e(line); p = line + 4;
        (void)strtod(p, &p);
        nav->crs = strtod(p, &p); nav->dn = strtod(p, &p); nav->M0 = strtod(p, &p);

        if (!fgets(line, sizeof(line), fp)) break; str_replace_d_e(line); p = line + 4;
        nav->cuc = strtod(p, &p); nav->e = strtod(p, &p); nav->cus = strtod(p, &p); nav->sqrtA = strtod(p, &p);

        if (!fgets(line, sizeof(line), fp)) break; str_replace_d_e(line); p = line + 4;
        nav->toe_sec = strtod(p, &p); nav->cic = strtod(p, &p); nav->OMEGA0 = strtod(p, &p); nav->cis = strtod(p, &p);

        if (!fgets(line, sizeof(line), fp)) break; str_replace_d_e(line); p = line + 4;
        nav->i0 = strtod(p, &p); nav->crc = strtod(p, &p); nav->omega = strtod(p, &p); nav->OMEGAdot = strtod(p, &p);

        if (!fgets(line, sizeof(line), fp)) break; str_replace_d_e(line); p = line + 4;
        nav->IDOT = strtod(p, &p);

        fgets(line, sizeof(line), fp); fgets(line, sizeof(line), fp);
        nav_count[idx]++;
        nav_total++;
    }
    fclose(fp);
    printf("已从 NAV 文件加载 %d 条星历。\n", nav_total);
}

void read_obs(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) { printf("错误: 无法打开 OBS 文件: %s\n", filename); return; }
    char line[4096];

    idx_G_D1 = idx_G_D2 = -1;
    idx_E_D1 = idx_E_D2 = -1;
    idx_C_D1 = idx_C_D2 = -1;

    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "APPROX POSITION XYZ")) {
            (void)sscanf(line, "%lf %lf %lf", &STA_POS_XYZ.x, &STA_POS_XYZ.y, &STA_POS_XYZ.z);
        }

        if (strstr(line, "SYS / # / OBS TYPES")) {
            char sys_code = line[0];
            if (sys_code != 'G' && sys_code != 'E' && sys_code != 'C') continue;

            char buffer[4096]; strcpy(buffer, line);
            long pos = ftell(fp); char next[1024];
            while (fgets(next, sizeof(next), fp)) {
                if (strstr(next, "SYS / # / OBS TYPES")) {
                    if (next[0] != sys_code && next[0] != ' ') { fseek(fp, pos, SEEK_SET); break; }
                    strcat(buffer, next); pos = ftell(fp);
                }
                else { fseek(fp, pos, SEEK_SET); break; }
            }

            char* token = strtok(buffer, " \t\n\r");
            int raw_count = 0, type_idx = 0;
            int* p_d1 = NULL;
            int* p_d2 = NULL;

            if (sys_code == 'G') { p_d1 = &idx_G_D1; p_d2 = &idx_G_D2; }
            else if (sys_code == 'E') { p_d1 = &idx_E_D1; p_d2 = &idx_E_D2; }
            else if (sys_code == 'C') { p_d1 = &idx_C_D1; p_d2 = &idx_C_D2; }

            while (token != NULL) {
                if (strcmp(token, "SYS") != 0 && strcmp(token, "/") != 0 && strcmp(token, "#") != 0 && strcmp(token, "OBS") != 0 && strcmp(token, "TYPES") != 0) {
                    if (raw_count >= 2) {
                        if (sys_code == 'G') {
                            if (strstr(token, "D1C")) *p_d1 = type_idx;
                            else if (strstr(token, "D1X") && *p_d1 == -1) *p_d1 = type_idx;
                            if (strstr(token, "D2W")) *p_d2 = type_idx;
                            else if (strstr(token, "D2P") && *p_d2 == -1) *p_d2 = type_idx;
                        }
                        else if (sys_code == 'E') {
                            if (strstr(token, "D1X")) *p_d1 = type_idx;
                            else if (strstr(token, "D1C") && *p_d1 == -1) *p_d1 = type_idx;
                            if (strstr(token, "D5X")) *p_d2 = type_idx;
                            else if (strstr(token, "D5Q") && *p_d2 == -1) *p_d2 = type_idx;
                        }
                        // [重点修正] 优先选择 D1I/D1X (B1C) 以匹配 1575.42 MHz
                        else if (sys_code == 'C') {
                            if (strstr(token, "D1I") || strstr(token, "D1X")) {
                                *p_d1 = type_idx;
                                cur_bds_f1 = FREQ_B1C; // 1575.42 MHz
                            }
                            else if (strstr(token, "D2I") && *p_d1 == -1) { // 仅当没有 D1I 时才退而求其次
                                *p_d1 = type_idx;
                                cur_bds_f1 = FREQ_B1I; // 1561.098 MHz
                            }

                            if (strstr(token, "D6I")) *p_d2 = type_idx; // B3I
                            else if (strstr(token, "D7I") && *p_d2 == -1) *p_d2 = type_idx;
                        }
                        type_idx++;
                    }
                    raw_count++;
                }
                token = strtok(NULL, " \t\n\r");
            }
        }
        if (strstr(line, "END OF HEADER")) break;
    }

    printf("索引: G(%d,%d) E(%d,%d) C(%d,%d)\n", idx_G_D1, idx_G_D2, idx_E_D1, idx_E_D2, idx_C_D1, idx_C_D2);
    if (idx_C_D1 != -1) printf("BDS D1 频率: %.3f MHz\n", cur_bds_f1 / 1e6);

    obs_count = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') {
            if (obs_count >= MAX_EPOCH) break;
            epoch_obs_t* ep = &obs_data[obs_count];
            double ep_time[6];
            int flag = 0, num_sats = 0;
            if (sscanf(line + 2, "%lf %lf %lf %lf %lf %lf %d %d",
                &ep_time[0], &ep_time[1], &ep_time[2], &ep_time[3], &ep_time[4], &ep_time[5],
                &flag, &num_sats) < 8) continue;

            ep->time = epoch2time(ep_time);

            if (flag != 0 && flag != 1) {
                for (int k = 0; k < num_sats; k++) fgets(line, sizeof(line), fp);
                continue;
            }

            int valid_sats = 0;
            for (int k = 0; k < num_sats; k++) {
                if (!fgets(line, sizeof(line), fp)) break;

                int sys = -1;
                if (line[0] == 'G') sys = SYS_GPS;
                else if (line[0] == 'E') sys = SYS_GAL;
                else if (line[0] == 'C') sys = SYS_BDS;
                else continue;

                int prn = atoi(line + 1);
                char field[17] = { 0 };
                double d1 = 0.0, d2 = 0.0;

                int i1 = -1, i2 = -1;
                if (sys == SYS_GPS) { i1 = idx_G_D1; i2 = idx_G_D2; }
                else if (sys == SYS_GAL) { i1 = idx_E_D1; i2 = idx_E_D2; }
                else if (sys == SYS_BDS) { i1 = idx_C_D1; i2 = idx_C_D2; }

                if (i1 >= 0 && strlen(line) > 3 + i1 * 16) { strncpy(field, line + 3 + i1 * 16, 14); d1 = strtod(field, NULL); }
                if (i2 >= 0 && strlen(line) > 3 + i2 * 16) { strncpy(field, line + 3 + i2 * 16, 14); memset(field, 0, 17); d2 = strtod(field, NULL); }

                if (d1 != 0.0) {
                    ep->sats[valid_sats].prn = prn;
                    ep->sats[valid_sats].sys = sys;
                    ep->sats[valid_sats].D1 = d1;
                    ep->sats[valid_sats].D2 = d2;
                    valid_sats++;
                }
            }
            ep->sat_num = valid_sats;
            if (valid_sats > 0) obs_count++;
        }
    }
    fclose(fp);
    printf("已加载 %d 个有效历元。\n", obs_count);
}

// 核心算法
void sat_pos_vel(gtime_t time, nav_t* nav, vec3_t* pos, vec3_t* vel, double* dts_rate) {
    double t_obs = time2gpst(time, NULL);

    // [修正1] 时间系统差异修正 (BDT = GPST - 14s)
    double time_bias = 0.0;
    if (nav->sys == SYS_BDS) {
        time_bias = -14.0;
    }

    double tk = (t_obs + time_bias) - nav->toe_sec;
    if (tk > 302400.0) tk -= 604800.0;
    if (tk < -302400.0) tk += 604800.0;

    double gm, omega_e_val;
    if (nav->sys == SYS_GAL) { gm = GM_GAL; omega_e_val = OMEGA_E; }
    else if (nav->sys == SYS_BDS) { gm = GM_BDS; omega_e_val = OMEGA_E_BDS; }
    else { gm = GM_GPS; omega_e_val = OMEGA_E; }

    double A = nav->sqrtA * nav->sqrtA;
    double n0 = sqrt(gm / (A * A * A));
    double n = n0 + nav->dn;
    double M = nav->M0 + n * tk;
    double E = M, E_old;
    for (int i = 0; i < 30; i++) { E_old = E; E = M + nav->e * sin(E); if (fabs(E - E_old) < 1e-13) break; }

    double sinE = sin(E), cosE = cos(E);
    double phi = atan2(sqrt(1.0 - nav->e * nav->e) * sinE, cosE - nav->e) + nav->omega;
    double u = phi + nav->cuc * cos(2.0 * phi) + nav->cus * sin(2.0 * phi);
    double r = A * (1.0 - nav->e * cosE) + nav->crc * cos(2.0 * phi) + nav->crs * sin(2.0 * phi);
    double i_orb = nav->i0 + nav->IDOT * tk + nav->cic * cos(2.0 * phi) + nav->cis * sin(2.0 * phi);

    // [修正2] GEO 卫星特殊处理 (PRN 1-5) 
    if (nav->sys == SYS_BDS && nav->prn <= 5) {
        // 1. 惯性系升交点经度
        double Omega_k = nav->OMEGA0 + nav->OMEGAdot * tk - omega_e_val * nav->toe_sec;

        // 2. 轨道面坐标
        double x_p = r * cos(u);
        double y_p = r * sin(u);

        // 3. 惯性系坐标 (GK)
        double cosO = cos(Omega_k), sinO = sin(Omega_k);
        double cosi = cos(i_orb), sini = sin(i_orb);
        double x_g = x_p * cosO - y_p * cosi * sinO;
        double y_g = x_p * sinO + y_p * cosi * cosO;
        double z_g = y_p * sini;

        // 4. 旋转到地固系: Rz(w*tk) * Rx(-5 deg) * X_g
        double alpha = -5.0 * PI / 180.0;
        double x_g2 = x_g;
        double y_g2 = y_g * cos(alpha) + z_g * sin(alpha);
        double z_g2 = -y_g * sin(alpha) + z_g * cos(alpha);

        double beta = omega_e_val * tk;
        pos->x = x_g2 * cos(beta) + y_g2 * sin(beta);
        pos->y = -x_g2 * sin(beta) + y_g2 * cos(beta);
        pos->z = z_g2;

        // 5. GEO 速度计算
        double E_dot = n / (1.0 - nav->e * cosE);
        double phi_dot = sqrt(1.0 - nav->e * nav->e) * E_dot / (1.0 - nav->e * cosE);
        double u_dot = phi_dot * (1.0 + 2.0 * (nav->cus * cos(2.0 * phi) - nav->cuc * sin(2.0 * phi)));
        double r_dot = A * nav->e * sinE * E_dot + 2.0 * (nav->crs * cos(2.0 * phi) - nav->crc * sin(2.0 * phi)) * phi_dot;
        double i_dot = nav->IDOT + 2.0 * (nav->cis * cos(2.0 * phi) - nav->cic * sin(2.0 * phi)) * phi_dot;
        double O_dot = nav->OMEGAdot;

        double x_p_dot = r_dot * cos(u) - r * u_dot * sin(u);
        double y_p_dot = r_dot * sin(u) + r * u_dot * cos(u);

        // dX_g/dt
        double dx_g = x_p_dot * cosO - x_p * sinO * O_dot - (y_p_dot * cosi * sinO + y_p * (-sini * i_dot * sinO + cosi * cosO * O_dot));
        double dy_g = x_p_dot * sinO + x_p * cosO * O_dot + (y_p_dot * cosi * cosO + y_p * (-sini * i_dot * cosO - cosi * sinO * O_dot));
        double dz_g = y_p_dot * sini + y_p * cosi * i_dot;

        // 旋转速度
        double rx_x = x_g;
        double rx_y = y_g * cos(alpha) + z_g * sin(alpha);

        double drx_x = dx_g;
        double drx_y = dy_g * cos(alpha) + dz_g * sin(alpha);
        double drx_z = -dy_g * sin(alpha) + dz_g * cos(alpha);

        double term1_x = drx_x * cos(beta) + drx_y * sin(beta);
        double term1_y = -drx_x * sin(beta) + drx_y * cos(beta);
        double term1_z = drx_z;

        double term2_x = omega_e_val * (-rx_x * sin(beta) + rx_y * cos(beta));
        double term2_y = omega_e_val * (-rx_x * cos(beta) - rx_y * sin(beta));

        vel->x = term1_x + term2_x;
        vel->y = term1_y + term2_y;
        vel->z = term1_z;

    }
    else {
        // --- 非 GEO 卫星 ---
        double L = nav->OMEGA0 + (nav->OMEGAdot - omega_e_val) * tk - omega_e_val * nav->toe_sec;
        double x_p = r * cos(u), y_p = r * sin(u);
        double cosL = cos(L), sinL = sin(L), cosi = cos(i_orb), sini = sin(i_orb);
        pos->x = x_p * cosL - y_p * cosi * sinL;
        pos->y = x_p * sinL + y_p * cosi * cosL;
        pos->z = y_p * sini;

        double E_dot = n / (1.0 - nav->e * cosE);
        double phi_dot = sqrt(1.0 - nav->e * nav->e) * E_dot / (1.0 - nav->e * cosE);
        double u_dot = phi_dot * (1.0 + 2.0 * (nav->cus * cos(2.0 * phi) - nav->cuc * sin(2.0 * phi)));
        double r_dot = A * nav->e * sinE * E_dot + 2.0 * (nav->crs * cos(2.0 * phi) - nav->crc * sin(2.0 * phi)) * phi_dot;
        double i_dot = nav->IDOT + 2.0 * (nav->cis * cos(2.0 * phi) - nav->cic * sin(2.0 * phi)) * phi_dot;
        double L_dot = nav->OMEGAdot - omega_e_val;

        double x_p_dot = r_dot * cos(u) - r * u_dot * sin(u);
        double y_p_dot = r_dot * sin(u) + r * u_dot * cos(u);

        vel->x = (x_p_dot * cosL - x_p * sinL * L_dot) - (y_p_dot * cosi * sinL + y_p * (-sini * i_dot * sinL + cosi * cosL * L_dot));
        vel->y = (x_p_dot * sinL + x_p * cosL * L_dot) + (y_p_dot * cosi * cosL + y_p * (-sini * i_dot * cosL - cosi * sinL * L_dot));
        vel->z = y_p_dot * sini + y_p * cosi * i_dot;
    }

    *dts_rate = nav->a1 + 2.0 * nav->a2 * tk;
}

int main() {
    printf("GNSS 多普勒测速程序 (GPS + Galileo + BDS)\n");
    read_nav("brdc1590.24p");
    read_obs("jfng1590.24o");
    FILE* fout = fopen("velocity_result.txt", "w");
    if (!fout) { printf("错误: 无法创建输出文件。\n"); return -1; }
    fprintf(fout, "Time(GPST)         Sys     Ns  PDOP   Ve(m/s)   Vn(m/s)   Vu(m/s)\n");

    int solved_epochs = 0;

    for (int i = 0; i < obs_count; i++) {
        epoch_obs_t* obs = &obs_data[i];

        int week_chk;
        double sec_chk = time2gpst(obs->time, &week_chk);
        double tod_chk = fmod(sec_chk, 86400.0);

        if (tod_chk < 7200.0 || tod_chk > 21600.0) continue;

        vec3_t rcv_vel = { 0 };
        double rcv_clk_drift = 0.0;
        int solved = 0;
        double final_pdop = 0.0;
        vec3_t final_enu = { 0 };
        int num_used = 0;

        for (int iter = 0; iter < 10; iter++) {
            double H[MAX_SAT * MAX_SYS * 4], y_vec[MAX_SAT * MAX_SYS], W[MAX_SAT * MAX_SYS];
            int n = 0;

            for (int j = 0; j < obs->sat_num; j++) {
                int prn = obs->sats[j].prn;
                int sys = obs->sats[j].sys;
                int idx = get_idx(prn, sys);

                if (nav_count[idx] == 0) continue;

                nav_t* best_nav = NULL;
                double min_dt = 1e9;
                for (int k = 0; k < nav_count[idx]; k++) {
                    double dt = fabs(timediff(obs->time, nav_data[idx][k].toc));
                    if (dt < min_dt) { min_dt = dt; best_nav = &nav_data[idx][k]; }
                }
                if (!best_nav || min_dt > 7200.0) continue;

                vec3_t sat_pos, sat_vel; double dts_rate;
                sat_pos_vel(obs->time, best_nav, &sat_pos, &sat_vel, &dts_rate);

                double r = sqrt(pow(sat_pos.x - STA_POS_XYZ.x, 2) + pow(sat_pos.y - STA_POS_XYZ.y, 2) + pow(sat_pos.z - STA_POS_XYZ.z, 2));
                double ex = (sat_pos.x - STA_POS_XYZ.x) / r, ey = (sat_pos.y - STA_POS_XYZ.y) / r, ez = (sat_pos.z - STA_POS_XYZ.z) / r;

                vec3_t diff = { sat_pos.x - STA_POS_XYZ.x, sat_pos.y - STA_POS_XYZ.y, sat_pos.z - STA_POS_XYZ.z };
                double E_enu, N_enu, U_enu;
                xyz2enu(STA_POS_XYZ, diff, &E_enu, &N_enu, &U_enu);
                if (asin(U_enu / r) < 10.0 * PI / 180.0) continue;

                double rr_obs = 0.0;
                double f1, f2, lam1, lam2;

                if (sys == SYS_GPS) { f1 = FREQ_L1; f2 = FREQ_L2; }
                else if (sys == SYS_GAL) { f1 = FREQ_E1; f2 = FREQ_E5a; }
                else {
                    // [重点] 使用自动识别的正确 BDS B1 频率 (1575.42 or 1561.098)
                    f1 = cur_bds_f1;
                    f2 = FREQ_B3I;
                }

                lam1 = CLIGHT / f1;
                lam2 = CLIGHT / f2;

                if (obs->sats[j].D1 != 0.0) {
                    if (obs->sats[j].D2 != 0.0) {
                        double rr1 = -obs->sats[j].D1 * lam1;
                        double rr2 = -obs->sats[j].D2 * lam2;
                        double gamma = (f1 * f1) / (f2 * f2);
                        rr_obs = (gamma * rr1 - rr2) / (gamma - 1.0);
                    }
                    else {
                        rr_obs = -obs->sats[j].D1 * lam1;
                    }
                }
                else continue;

                double v_sat_proj = ex * sat_vel.x + ey * sat_vel.y + ez * sat_vel.z;
                double v_rcv_proj = ex * rcv_vel.x + ey * rcv_vel.y + ez * rcv_vel.z;
                double sagnac_corr = (OMEGA_E / CLIGHT) * (sat_vel.x * STA_POS_XYZ.y - sat_vel.y * STA_POS_XYZ.x);
                double rr_comp = v_sat_proj - v_rcv_proj + CLIGHT * dts_rate - rcv_clk_drift + sagnac_corr;

                H[n * 4 + 0] = -ex; H[n * 4 + 1] = -ey; H[n * 4 + 2] = -ez; H[n * 4 + 3] = 1.0;
                y_vec[n] = rr_obs - rr_comp;

                double el = asin(U_enu / r);
                double el_deg = el * 180.0 / PI;
                double k_sys = 1.0;
                double weight_base = (el_deg >= 30.0) ? 1.0 : (4.0 * sin(el) * sin(el));
                if (el_deg <= 0.0) weight_base = 0.0001;

                W[n] = weight_base * k_sys * k_sys;
                n++;
            }

            if (n < 4) break;

            double AtWA[16] = { 0 }, AtWv[4] = { 0 };
            for (int k = 0; k < n; k++) {
                for (int r = 0; r < 4; r++) {
                    for (int c = 0; c < 4; c++) AtWA[r * 4 + c] += H[k * 4 + r] * W[k] * H[k * 4 + c];
                    AtWv[r] += H[k * 4 + r] * W[k] * y_vec[k];
                }
            }
            double inv[16], dx[4] = { 0 };
            if (!matinv4(AtWA, inv)) break;
            for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) dx[r] += inv[r * 4 + c] * AtWv[c];

            rcv_vel.x += dx[0]; rcv_vel.y += dx[1]; rcv_vel.z += dx[2]; rcv_clk_drift += dx[3];

            if (sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) < 1e-4) {
                final_pdop = sqrt(inv[0] + inv[5] + inv[10]);
                vec3_t vres = { rcv_vel.x, rcv_vel.y, rcv_vel.z };
                xyz2enu(STA_POS_XYZ, vres, &final_enu.x, &final_enu.y, &final_enu.z);
                solved = 1; num_used = n;
                break;
            }
        }

        if (solved) {
            int h = (int)(tod_chk / 3600), m = (int)((tod_chk - h * 3600) / 60);
            double s = tod_chk - h * 3600 - m * 60;
            fprintf(fout, "%02d:%02d:%02.0f       GPS/GAL/BDS %2d %5.2f %8.3f  %8.3f  %8.3f\n",
                h, m, s, num_used, final_pdop, final_enu.x, final_enu.y, final_enu.z);
            solved_epochs++;
        }
    }
    printf("处理完成。共解算 %d 个历元。\n", solved_epochs);
    if (solved_epochs == 0) printf("错误: 没有解算任何历元。\n");
    fclose(fout);
    return 0;
}