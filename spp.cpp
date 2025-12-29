#define _CRT_SECURE_NO_WARNINGS  // 防止 Visual Studio 报 strcpy 等不安全函数的警告
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
using namespace std;
// GPS时间结构体：周 + 周内秒
struct GPSTime {
    int week;  // GPS周
    double sec;// 周内秒 (0 - 604800)
};
// 导航电文数据结构 存储广播星历
struct NavData {
    char sys;       // 卫星系统 (G/C/E/R)
    int sat;        // 卫星号
    GPSTime toe;    // 星历参考时间(Time of Ephemeris)
    GPSTime toc;    // 钟差参考时间(Time of Clock)

    // 开普勒轨道参数
    double M0;      // 参考时刻的平近点角
    double e;       // 轨道偏心率
    double sqrtA;   // 长半轴的平方根
    double Omega0;  // 参考时刻升交点赤经
    double i0;      // 参考时刻轨道倾角
    double omega;   // 近地点角距
    double cuc, cus;// 升交点角距的余弦/正弦调和改正项
    double crc, crs;// 轨道半径的余弦/正弦调和改正项
    double cic, cis;// 轨道倾角的余弦/正弦调和改正项
    double OmegaDot;// 升交点赤经变化率
    double IDot;    // 轨道倾角变化率
    double delta_n; // 平均角速度的修正量

    // 钟差参数
    double a0, a1, a2;//钟差钟速钟漂
    double tgd;     // 群延迟 (Total Group Delay)，用于修正不同频率间的硬件延迟

    // GLONASS 特定，后续增加
    double X, Y, Z;      // 位置 (PZ-90坐标系)
    double Vx, Vy, Vz;  // Vel
    double Ax, Ay, Az;  /// 加速度 (主要受日月引力摄动)
    int freq_num;       // 频率号 (GLONASS是频分多址 FDMA)
};
// 观测数据结构 (存储.O文件读入的数据)
struct ObsData {
    GPSTime time; // 观测时刻 (接收机时间)
    char sys;
    int sat;
    double P1, P2; //双频伪距观测值 (单位: 米)
};

const double PI = 3.1415926535897932;
const double CLIGHT = 299792458.0;
// 地球引力常数 (GM)
const double GM_GPS = 3.9860050E14;
const double GM_BDS = 3.986004418E14;
const double GM_GAL = 3.986004418E14;
const double GM_GLO = 3.9860044E14;
// 地球自转角速度 (用于sagnac效应修正和坐标系旋转)
const double OMEGA_GPS = 7.2921151467E-5;
const double OMEGA_BDS = 7.2921150E-5;
const double OMEGA_GAL = 7.2921151467E-5;
const double OMEGA_GLO = 7.2921150E-5;
const double OMEGA_E = 7.2921151467E-5; // 通用的地球自转角速度  
// WGS84 椭球参数
const double RE_WGS84 = 6378137.0;     //长半轴
const double FE_WGS84 = 1.0 / 298.257223563; // 扁率
const double LEAPS = 18.0; // 当前跳秒(GPS与UTC的差)

// 载波频率定义 (用于无电离层组合计算)
const double FREQ_GPS_L1 = 1.57542E9;
const double FREQ_GPS_L2 = 1.22760E9;
const double FREQ_BDS_B1 = 1.561098E9;
const double FREQ_BDS_B3 = 1.26852E9;// 北斗通常用 B1/B3 组合
const double FREQ_GAL_E1 = 1.57542E9;
const double FREQ_GAL_E5a = 1.17645E9;
const double FREQ_GLO_G1_BASE = 1.60200E9;
const double FREQ_GLO_G1_DELTA = 0.5625E6;
const double FREQ_GLO_G2_BASE = 1.24600E9;
const double FREQ_GLO_G2_DELTA = 0.4375E6;

// SPP 截止高度角 (低于15度的卫星不参与解算)
const double EL_MASK = 15.0 * PI / 180.0;

// 函数原型声明！
// gtime
GPSTime epoch2time(const double* ep);  //把年、月、日、时、分、秒（ epoch）转换成 GPS 时间（周和周内秒 ）。
void time2epoch(GPSTime t, double* ep); //把 GPS 时间（周和周内秒 ）转换成 年、月、日、时、分、秒（ epoch）。
double timediff(GPSTime t1, GPSTime t2); //计算两个 GPS 时间点之间的时间差
GPSTime gpst2utc(GPSTime t);     //GPS 时间和 UTC 之间进行转换
GPSTime utc2gpst(GPSTime t);

// matrix
void matmul(const char* tr, int n, int k, int m, double alpha, const double* A, const double* B, double beta, double* C); //矩阵乘法
int matinv(double* A, int n);  //矩阵求逆

// tropo
int trop_model_prec(GPSTime time, const double* pos, const double* azel, double* trop, double* var);//对流层延迟

// orbit
void satPos(GPSTime tTx, const NavData& nav, double* rs, double* dts);//根据导航电文，计算卫星在信号发射时刻的位置 (x, y, z) 和卫星钟差 dts。
void satPosVel(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift); //多一个速度

// rinex
std::vector<NavData> readNavFile(const std::string& filename);  //读.p
std::vector<std::vector<ObsData>> readObsFile(const std::string& filename, double* approxPos);//读.o

// spp 核心
void spp_gnss_process(const std::vector<std::vector<ObsData>>& obsList, const std::vector<NavData>& navList, const double* refPos, const std::string& outfile);

// Common Utils
void ecef2pos_common(const double* r, double* pos);
void xyz2enu_common(const double* r, const double* ref, double* enu);
void sat_azel_common(const double* r, const double* rs, double* azel);//计算卫星相对于接收机的方位角和高度角



// 模块实现: 时间处理 
static bool isLeapYear(int year) {   //判断闰年
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

static int getDayOfYear(int year, int month, int day) {  //算是一年中第几天
    int days_in_month[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    if (isLeapYear(year)) days_in_month[1] = 29;
    int doy = 0;
    for (int i = 0; i < month - 1; i++) doy += days_in_month[i];
    doy += day;
    return doy;
}

static int daysBetween(int year1, int month1, int day1, int year2, int month2, int day2) { //计算天数差
    long long days1 = (long long)year1 * 365 + getDayOfYear(year1, month1, day1);
    long long days2 = (long long)year2 * 365 + getDayOfYear(year2, month2, day2);
    int y_start = std::min(year1, year2);
    int y_end = std::max(year1, year2);
    for (int y = y_start; y < y_end; y++) {
        if (isLeapYear(y)) {
            if (y < year2) days2++;
            else days1++;
        }
    }
    return (int)(days2 - days1);
}

GPSTime epoch2time(const double* ep) {
    GPSTime t = { 0, 0.0 };
    int year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];
    int hour = (int)ep[3], min = (int)ep[4];
    double sec = ep[5];
    int days_from_gps = daysBetween(1980, 1, 6, year, mon, day);
    t.week = days_from_gps / 7;
    int day_of_week = days_from_gps % 7;
    if (day_of_week < 0) { day_of_week += 7; t.week -= 1; }
    t.sec = day_of_week * 86400.0 + hour * 3600.0 + min * 60.0 + sec;
    return t;
}

void time2epoch(GPSTime t, double* ep) {
    double jd = t.week * 7.0 + t.sec / 86400.0 + 2444244.5;
    double z = floor(jd + 0.5);
    double f = jd + 0.5 - z;
    double alpha = floor((z - 1867216.25) / 36524.25);
    double a = z + 1 + alpha - floor(alpha / 4.0);
    double b = a + 1524;
    double c = floor((b - 122.1) / 365.25);
    double d = floor(365.25 * c);
    double e = floor((b - d) / 30.6001);
    ep[2] = b - d - floor(30.6001 * e) + f;
    ep[1] = e < 14 ? e - 1 : e - 13;
    ep[0] = ep[1] > 2 ? c - 4716 : c - 4715;
    double day_int = floor(ep[2]);
    double day_frac = ep[2] - day_int;
    ep[2] = day_int;
    double total_sec = day_frac * 86400.0;
    ep[3] = floor(total_sec / 3600.0);
    total_sec -= ep[3] * 3600.0;
    ep[4] = floor(total_sec / 60.0);
    ep[5] = total_sec - ep[4] * 60.0;
    if (ep[5] >= 60.0) { ep[5] = 0.0; ep[4] += 1; }
    if (ep[4] >= 60.0) { ep[4] = 0.0; ep[3] += 1; }
}

double timediff(GPSTime t1, GPSTime t2) {
    return (t1.week - t2.week) * 604800.0 + (t1.sec - t2.sec);
}

static double str2num_char(const char* s, int i, int n) {
    char str[256], * p = str;
    if (i < 0 || n < 1 || sizeof(str) < n + 1) return 0.0;
    for (s += i; *s && --n >= 0; s++) *p++ = *s == 'D' || *s == 'd' ? 'E' : *s;
    *p = '\0';
    return atof(str);
}

//模块实现: 矩阵运算
// 矩阵乘法: C = alpha * A * B + beta * C
// tr: 字符串，"NN" 表示 A不转置B不转置，"TN" 表示 A转置B不转置 
void matmul(const char* tr, int n, int k, int m, double alpha,
    const double* A, const double* B, double beta, double* C)
{
    double d;
    int f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            d = 0.0;
            for (int x = 0; x < k; x++) {
                double a = (f == 1 || f == 2) ? A[i * k + x] : A[x * n + i];
                double b = (f == 1 || f == 3) ? B[x * m + j] : B[j * k + x];
                d += a * b;
            }
            if (beta == 0.0) C[i * m + j] = alpha * d;
            else C[i * m + j] = alpha * d + beta * C[i * m + j];
        }
    }
}
// 矩阵求逆 (高斯-约旦消元法)
// 输入 A，输出 A 的逆矩阵覆盖原数组
int matinv(double* A, int n) {
    int i, j, k;
    vector<int> indxc(n), indxr(n), ipiv(n, 0);
    for (i = 0; i < n; i++) {
        double big = 0.0;
        int irow = -1, icol = -1;
        for (j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(A[j * n + k]) >= big) {
                            big = fabs(A[j * n + k]); irow = j; icol = k;
                        }
                    }
                }
            }
        }
        if (irow == -1 || icol == -1 || big == 0.0) return -1;
        ipiv[icol]++;
        if (irow != icol) {
            for (int l = 0; l < n; l++) swap(A[irow * n + l], A[icol * n + l]);
        }
        indxr[i] = irow; indxc[i] = icol;
        if (A[icol * n + icol] == 0.0) return -1;
        double pivinv = 1.0 / A[icol * n + icol];
        A[icol * n + icol] = 1.0;
        for (int l = 0; l < n; l++) A[icol * n + l] *= pivinv;
        for (int ll = 0; ll < n; ll++) {
            if (ll != icol) {
                double dum = A[ll * n + icol];
                A[ll * n + icol] = 0.0;
                for (int l = 0; l < n; l++) A[ll * n + l] -= A[icol * n + l] * dum;
            }
        }
    }
    for (int l = n - 1; l >= 0; l--) {
        if (indxr[l] != indxc[l]) {
            for (int k = 0; k < n; k++) swap(A[k * n + indxr[l]], A[k * n + indxc[l]]);
        }
    }
    return 0;
}

// 模块实现: 对流层模型 

// 计算年积日 (1.0 - 366.0)，用于对流层季节性参数变化
static double time2doy_trop(GPSTime t) {
    double ep[6]; time2epoch(t, ep);
    double epoch_first[6] = { ep[0], 1, 1, 0, 0, 0 };
    GPSTime t_first = epoch2time(epoch_first);
    double dt = timediff(t, t_first);
    return dt / 86400.0 + 1.0;
}

// 标准气象模型: 根据高度估算气压(P)、温度(T)、水汽压(e)
static void get_met_standard(double h, double* pres, double* temp, double* e) {
    const double PR_STD = 1013.25;
    const double TR_STD = 288.15;
    const double HR_STD = 50.0;
    if (h < -100.0) h = -100.0;
    if (h > 10000.0) h = 10000.0;
    *pres = PR_STD * pow(1.0 - 0.0000226 * h, 5.225);
    *temp = TR_STD - 0.0065 * h;
    double humidity = HR_STD / 100.0;
    double Tc = *temp - 273.15;
    double es = 6.11 * pow(10.0, (7.5 * Tc) / (237.3 + Tc));
    *e = humidity * es;
}

// Hopfield 模型: 计算天顶方向的干延迟(zhd)和湿延迟(zwd)
static void trop_hopfield_zenith(double P, double T, double e, double* zhd, double* zwd) {
    double h_dry = 40136.0 + 148.72 * (T - 273.16);  //干分量等效高度
    double h_wet = 11000.0;                         // 湿分量等效高度
    double N_d0 = 77.64 * (P / T);                  // 地面干折射率
    double N_w0 = -12.96 * (e / T) + 3.718e5 * (e / (T * T)); //地面湿折射率
    *zhd = 1.0e-6 / 5.0 * N_d0 * h_dry;
    *zwd = 1.0e-6 / 5.0 * N_w0 * h_wet;
}
// NMF映射函数的辅助公式
static double nmf_func(double el, double a, double b, double c) {
    double sin_el = sin(el);
    double top = 1.0 + a / (1.0 + b / (1.0 + c));
    double bot = sin_el + a / (sin_el + b / (sin_el + c));
    return top / bot;
}
// NMF 系数表 (按纬度分布)，查表用
static const double COEF_DRY[5][3][2] = {
    {{1.2769934e-3, 0.0}, {2.9153695e-3, 0.0}, {62.610505e-3, 0.0}},
    {{1.2683230e-3, 1.2709626e-5}, {2.9152299e-3, 2.1414979e-5}, {62.837393e-3, 9.0128400e-5}},
    {{1.2465397e-3, 2.6523662e-5}, {2.9288445e-3, 3.0160779e-5}, {63.721774e-3, 4.3497037e-5}},
    {{1.2196049e-3, 3.4000452e-5}, {2.9022565e-3, 7.2543778e-5}, {63.824265e-3, 8.4795348e-5}},
    {{1.2045996e-3, 4.1202191e-5}, {2.9024912e-3, 11.723375e-5}, {64.258455e-3, 17.037206e-5}}
};
static const double COEF_WET[5][3] = {
    {5.8021897e-4, 1.4275268e-3, 4.3472964e-2},
    {5.6794847e-4, 1.5138625e-3, 4.6729510e-2},
    {5.8118019e-4, 1.4572752e-3, 4.3908931e-2},
    {5.9727542e-4, 1.5007428e-3, 4.4626982e-2},
    {6.1641693e-4, 1.7599082e-3, 5.4736038e-2}
};
// 计算 NMF 映射函数 ，将天顶延迟映射到当前卫星高度角路径
static void trop_nmf(GPSTime time, const double* pos, double el, double* map_dry, double* map_wet) {
    double lat = fabs(pos[0] * 180.0 / PI);
    double h = pos[2];
    double doy = time2doy_trop(time);
    double coef_d[3], coef_w[3];
    double lats[] = { 15.0, 30.0, 45.0, 60.0, 75.0 };

    // 干分量
    for (int k = 0; k < 3; k++) {
        double val_avg, val_amp;
        if (lat <= 15.0) { val_avg = COEF_DRY[0][k][0]; val_amp = COEF_DRY[0][k][1]; }
        else if (lat >= 75.0) { val_avg = COEF_DRY[4][k][0]; val_amp = COEF_DRY[4][k][1]; }
        else {
            int idx = (int)((lat - 15.0) / 15.0);
            double ratio = (lat - lats[idx]) / 15.0;
            val_avg = COEF_DRY[idx][k][0] + (COEF_DRY[idx + 1][k][0] - COEF_DRY[idx][k][0]) * ratio;
            val_amp = COEF_DRY[idx][k][1] + (COEF_DRY[idx + 1][k][1] - COEF_DRY[idx][k][1]) * ratio;
        }
        double sign = (pos[0] >= 0) ? 1.0 : -1.0;
        coef_d[k] = val_avg - sign * val_amp * cos(2.0 * PI * (doy - 28.0) / 365.25);
    }
    // 湿分量
    for (int k = 0; k < 3; k++) {
        if (lat <= 15.0) coef_w[k] = COEF_WET[0][k];
        else if (lat >= 75.0) coef_w[k] = COEF_WET[4][k];
        else {
            int idx = (int)((lat - 15.0) / 15.0);
            double ratio = (lat - lats[idx]) / 15.0;
            coef_w[k] = COEF_WET[idx][k] + (COEF_WET[idx + 1][k] - COEF_WET[idx][k]) * ratio;
        }
    }

    const double a_ht = 2.53e-5, b_ht = 5.49e-3, c_ht = 1.14e-3;
    double m_dry_0 = nmf_func(el, coef_d[0], coef_d[1], coef_d[2]);
    double m_ht = nmf_func(el, a_ht, b_ht, c_ht);
    *map_dry = m_dry_0 + (1.0 / sin(el) - m_ht) * (h / 1000.0);
    *map_wet = nmf_func(el, coef_w[0], coef_w[1], coef_w[2]);
}
// 对流层主函数
// 输入: 时间, 接收机位置, 卫星方位/高度角
// 输出: 总对流层延迟trop
int trop_model_prec(GPSTime time, const double* pos, const double* azel, double* trop, double* var) {
    double el = azel[1];
    if (el < 3.0 * PI / 180.0) { *trop = 0.0; if (var) *var = 0.0; return 0; }
    double pres, temp, e;
    get_met_standard(pos[2], &pres, &temp, &e);
    double zhd, zwd;
    trop_hopfield_zenith(pres, temp, e, &zhd, &zwd);
    double map_dry, map_wet;
    trop_nmf(time, pos, el, &map_dry, &map_wet);
    *trop = zhd * map_dry + zwd * map_wet;
    if (var) *var = 0.3 * 0.3;
    return 1;
}

// 模块实现: 轨道计算 

const double J2_GLO = 1.0826257E-3;
const double RE_GLO = 6378136.0;
const double SIN_5 = -0.0871557427476582;// -5度正弦 (北斗GEO旋转用)
const double COS_5 = 0.9961946980917456;
// 开普勒轨道计算 
static void satPosKepler(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift) {
    // 确定常数
    double GM_VAL, OMEGA_VAL;
    if (nav.sys == 'C') { GM_VAL = GM_BDS; OMEGA_VAL = OMEGA_BDS; }
    else if (nav.sys == 'E') { GM_VAL = GM_GAL; OMEGA_VAL = OMEGA_GAL; }
    else { GM_VAL = GM_GPS; OMEGA_VAL = OMEGA_GPS; }
    // 计算归化时间 tk
    double tk = timediff(tTx, nav.toe);
    if (nav.sys == 'C') tk -= 14.0; // BDS时间与GPS时间的14秒差
    // 计算平近点角 M
    double A = nav.sqrtA * nav.sqrtA;
    double n0 = sqrt(GM_VAL / (A * A * A));
    double n = n0 + nav.delta_n;
    double M = nav.M0 + n * tk;
    // 迭代解开普勒方程计算偏近点角 E: M = E - e*sin(E)
    double E = M, Eold;
    for (int i = 0; i < 30; i++) {
        Eold = E;
        E = M + nav.e * sin(E);
        if (fabs(E - Eold) < 1e-13) break;
    }
    double sinE = sin(E), cosE = cos(E);
    double Edot = n / (1.0 - nav.e * cosE);
    // 计算真近点角 phi
    double phi = atan2(sqrt(1.0 - nav.e * nav.e) * sinE, cosE - nav.e) + nav.omega;
    double sin2phi = sin(2.0 * phi), cos2phi = cos(2.0 * phi);
    //计算摄动改正项 (二阶调和)
    double du = nav.cus * sin2phi + nav.cuc * cos2phi;// 纬度幅角改正
    double dr = nav.crs * sin2phi + nav.crc * cos2phi;// 径向距离改正
    double di = nav.cis * sin2phi + nav.cic * cos2phi;// 轨道倾角改正
    // 计算修正后的参数
    double u = phi + du;
    double r = A * (1.0 - nav.e * cosE) + dr;
    double i = nav.i0 + nav.IDot * tk + di;
    // 计算轨道平面坐标
    double x_prime = r * cos(u);
    double y_prime = r * sin(u);
    //速度计算逻位置求导
    double nudot = Edot * sqrt(1.0 - nav.e * nav.e) / (1.0 - nav.e * cosE);
    double udot = nudot + 2.0 * nudot * (nav.cus * cos2phi - nav.cuc * sin2phi);
    double rdot = A * nav.e * sinE * Edot + 2.0 * nudot * (nav.crs * cos2phi - nav.crc * sin2phi);
    double idot = nav.IDot + 2.0 * nudot * (nav.cis * cos2phi - nav.cic * sin2phi);
    double vx_prime = rdot * cos(u) - r * sin(u) * udot;
    double vy_prime = rdot * sin(u) + r * cos(u) * udot;
    // 计算轨道平面坐标
    bool is_bds_geo = (nav.sys == 'C' && (nav.sat <= 5 || nav.sat >= 59));
    double x_final, y_final, z_final;
    double vx_final, vy_final, vz_final;
    double Omega, dOmega;

    if (is_bds_geo) {// 北斗 GEO 卫星特殊处理 (增加了5度倾角的旋转)
        Omega = nav.Omega0 + nav.OmegaDot * tk - OMEGA_VAL * nav.toe.sec;
        dOmega = nav.OmegaDot;
        double sinO = sin(Omega), cosO = cos(Omega);
        double sini = sin(i), cosi = cos(i);
        // 先计算惯性系坐标
        double xg = x_prime * cosO - y_prime * cosi * sinO;
        double yg = x_prime * sinO + y_prime * cosi * cosO;
        double zg = y_prime * sini;
       
        double vxg = vx_prime * cosO - vy_prime * cosi * sinO - (x_prime * sinO + y_prime * cosi * cosO) * dOmega + y_prime * sinO * sini * idot;
        double vyg = vx_prime * sinO + vy_prime * cosi * cosO + (x_prime * cosO - y_prime * cosi * sinO) * dOmega - y_prime * cosO * sini * idot;
        double vzg = vy_prime * sini + y_prime * cosi * idot;
        // 执行 RX(-5 deg) 和 RZ(wt) 旋转
        double rx = xg;
        double ry = yg * COS_5 + zg * SIN_5;
        double rz = -yg * SIN_5 + zg * COS_5;
        double alpha = OMEGA_VAL * tk;
        double sinA = sin(alpha), cosA = cos(alpha);

        x_final = rx * cosA + ry * sinA;
        y_final = -rx * sinA + ry * cosA;
        z_final = rz;

        double vrx = vxg;
        double vry = vyg * COS_5 + vzg * SIN_5;
        double vrz = -vyg * SIN_5 + vzg * COS_5;
        double dPx = rx * (-sinA * OMEGA_VAL) + ry * (cosA * OMEGA_VAL);
        double dPy = -rx * (cosA * OMEGA_VAL) + ry * (-sinA * OMEGA_VAL);

        vx_final = dPx + (vrx * cosA + vry * sinA);
        vy_final = dPy + (-vrx * sinA + vry * cosA);
        vz_final = vrz;

    }
    else {    
        // Omega 计算需减去地球自转角 (OMEGA_E * tk)，实现瞬间惯性系到地固系的转换
        Omega = nav.Omega0 + (nav.OmegaDot - OMEGA_VAL) * tk - OMEGA_VAL * nav.toe.sec;
        dOmega = nav.OmegaDot - OMEGA_VAL;
        double sinO = sin(Omega), cosO = cos(Omega);
        double sini = sin(i), cosi = cos(i);

        x_final = x_prime * cosO - y_prime * cosi * sinO;
        y_final = x_prime * sinO + y_prime * cosi * cosO;
        z_final = y_prime * sini;

        vx_final = (vx_prime * cosO - vy_prime * cosi * sinO) - (x_prime * sinO + y_prime * cosi * cosO) * dOmega + y_prime * sinO * sini * idot;
        vy_final = (vx_prime * sinO + vy_prime * cosi * cosO) + (x_prime * cosO - y_prime * cosi * sinO) * dOmega - y_prime * cosO * sini * idot;
        vz_final = vy_prime * sini + y_prime * cosi * idot;
    }
    //输出结果
    rs[0] = x_final; rs[1] = y_final; rs[2] = z_final;
    if (vs) { vs[0] = vx_final; vs[1] = vy_final; vs[2] = vz_final; }
    //卫星钟差计算  包含相对论效应
    if (dts) {
        *dts = nav.a0 + nav.a1 * tk + nav.a2 * tk * tk;
        // 相对论效应修正 (轨道偏心率引起的周期性变化)
        *dts -= 2.0 * sqrt(GM_VAL * A) * nav.e * sinE / (CLIGHT * CLIGHT);
    }
    if (dts_drift) *dts_drift = nav.a1 + 2.0 * nav.a2 * tk;
}
// GLONASS 微分方程 (位置和速度的导数)
static void glonass_deq(const double* x, double* xdot, const double* acc) {
    // 计算 r^2, r^3 等
    double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    double r3 = r2 * sqrt(r2), r5 = r2 * r2 * sqrt(r2);
    // J2 项摄动 (地球扁率影响)
    double a = 1.5 * J2_GLO * GM_GLO * RE_GLO * RE_GLO / r5;
    double b = 5.0 * x[2] * x[2] / r2;
    double c = -GM_GLO / r3 - a * (1.0 - b);
    // 速度导数 = 加速度
    xdot[0] = x[3]; xdot[1] = x[4]; xdot[2] = x[5];
    double omg2 = OMEGA_GLO * OMEGA_GLO;
    // 运动方程 (含地球自转离心力和科里奥利力项)
    xdot[3] = (c + omg2) * x[0] + 2.0 * OMEGA_GLO * x[4] + acc[0];
    xdot[4] = (c + omg2) * x[1] - 2.0 * OMEGA_GLO * x[3] + acc[1];
    xdot[5] = (c - 2.0 * a) * x[2] + acc[2];
}
// GLONASS 轨道积分 (使用 4阶 Runge-Kutta 方法)
// 因为GLONASS广播星历给出的是某一时刻的状态向量(P,V,A)，求其他时刻需要数值积分
static void glonass_orbit(double t, double* x, const double* acc) {
    double k1[6], k2[6], k3[6], k4[6], w[6];
    double step = (t < 0.0) ? -30.0 : 30.0;// 积分步长 30秒
    double dt = t;
    while (fabs(dt) > 1e-9) {    // 分步积分直到目标时间
        if (fabs(dt) < fabs(step)) step = dt;
        glonass_deq(x, k1, acc); for (int i = 0; i < 6; i++) w[i] = x[i] + k1[i] * step / 2.0;
        glonass_deq(w, k2, acc); for (int i = 0; i < 6; i++) w[i] = x[i] + k2[i] * step / 2.0;
        glonass_deq(w, k3, acc); for (int i = 0; i < 6; i++) w[i] = x[i] + k3[i] * step;
        glonass_deq(w, k4, acc); for (int i = 0; i < 6; i++) x[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * step / 6.0;
        dt -= step;
    }
}
// 统一的卫星位置计算入口
void satPosVel(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift) {
    if (nav.sys == 'R') {
        double x[6] = { nav.M0, nav.e, nav.sqrtA, nav.Omega0, nav.i0, nav.omega };
        double acc[3] = { nav.cuc, nav.cus, nav.crc };
        double t = timediff(tTx, nav.toe);
        glonass_orbit(t, x, acc);
        if (rs) { rs[0] = x[0]; rs[1] = x[1]; rs[2] = x[2]; }
        if (vs) { vs[0] = x[3]; vs[1] = x[4]; vs[2] = x[5]; }
        if (dts) *dts = -nav.a0 + nav.a1 * t;
        if (dts_drift) *dts_drift = nav.a1;
    }
    else {
        satPosKepler(tTx, nav, rs, vs, dts, dts_drift);
    }
}

void satPos(GPSTime tTx, const NavData& nav, double* rs, double* dts) {
    satPosVel(tTx, nav, rs, nullptr, dts, nullptr);
}

// 模块实现: RINEX 读取 
// 
// 辅助结构: 记录 RINEX 文件中不同信号所在的列号
// 例如: GPS的 C1C 在第0列, C2W 在第5列
struct ColIdx {
    int p1 = -1, p2 = -1, d1 = -1, d2 = -1;
    std::string c1, c2;
};
// 信号优先级表: 决定优先读取哪些信号进行解算 (如优先用 C1C/C2W)
static std::map<char, std::vector<std::string>> prio_P1 = {
    {'G', {"C1C", "C1W", "C1X"}}, {'E', {"C1X", "C1C"}},
    {'C', {"C2I", "C1P", "C1X"}}, {'R', {"C1C", "C1P"}}// 注意: 北斗B1通常对应 RINEX 3.04 中的 C2I (1561MHz)
};
static std::map<char, std::vector<std::string>> prio_P2 = {
    {'G', {"C2W", "C2P", "C2X"}}, {'E', {"C5X", "C5Q"}},
    {'C', {"C6I", "C7I"}},        {'R', {"C2C", "C2P"}}// 北斗B3通常对应 C6I
};

static double str2num_str(const std::string& s) {//D和E
    try {
        if (s.find_first_not_of(" \t\r\n") == std::string::npos) return 0.0;
        std::string temp = s;
        for (char& c : temp) if (c == 'D' || c == 'd') c = 'E';
        return stod(temp);
    }
    catch (...) { return 0.0; }
}

// 读取导航电文文件
std::vector<NavData> readNavFile(const std::string& filename) {
    std::vector<NavData> navList;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法打开导航文件 " << filename << std::endl;
        return navList;
    }
    std::string line;
    while (getline(file, line)) if (line.find("END OF HEADER") != std::string::npos) break;

    NavData nav = { 0 };
    std::vector<std::string> buff;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] != ' ' && line.length() > 5) {
            if (nav.sat != 0) navList.push_back(nav);
            nav = { 0 }; nav.sys = line[0];
            try { nav.sat = stoi(line.substr(1, 2)); }
            catch (...) { continue; }
            try {
                int yr = stoi(line.substr(4, 4));
                int mon = stoi(line.substr(9, 2));
                int day = stoi(line.substr(12, 2));
                int hr = stoi(line.substr(15, 2));
                int min = stoi(line.substr(18, 2));
                double sec = str2num_str(line.substr(21, 4));
                double toc_ep[6] = { (double)yr, (double)mon, (double)day, (double)hr, (double)min, sec };
                nav.toc = epoch2time(toc_ep);
            }
            catch (...) { continue; }
            nav.a0 = str2num_str(line.substr(23, 19));
            nav.a1 = str2num_str(line.substr(42, 19));
            nav.a2 = str2num_str(line.substr(61, 19));
            buff.clear();
        }
        else {
            buff.push_back(line);
            if (nav.sys == 'R') {
                if (buff.size() >= 3) {
                    nav.M0 = str2num_str(buff[0].substr(4, 19)) * 1000.0;
                    nav.Omega0 = str2num_str(buff[0].substr(23, 19)) * 1000.0;
                    nav.cuc = str2num_str(buff[0].substr(42, 19)) * 1000.0;
                    nav.e = str2num_str(buff[1].substr(4, 19)) * 1000.0;
                    nav.i0 = str2num_str(buff[1].substr(23, 19)) * 1000.0;
                    nav.cus = str2num_str(buff[1].substr(42, 19)) * 1000.0;
                    if (buff[1].length() >= 61) nav.freq_num = (int)str2num_str(buff[1].substr(61, 19));
                    nav.sqrtA = str2num_str(buff[2].substr(4, 19)) * 1000.0;
                    nav.omega = str2num_str(buff[2].substr(23, 19)) * 1000.0;
                    nav.crc = str2num_str(buff[2].substr(42, 19)) * 1000.0;
                    nav.toe = nav.toc;
                }
            }
            else {
                if (buff.size() >= 7) {
                    nav.crs = str2num_str(buff[0].substr(23, 19));
                    nav.delta_n = str2num_str(buff[0].substr(42, 19));
                    nav.M0 = str2num_str(buff[0].substr(61, 19));
                    nav.cuc = str2num_str(buff[1].substr(4, 19));
                    nav.e = str2num_str(buff[1].substr(23, 19));
                    nav.cus = str2num_str(buff[1].substr(42, 19));
                    nav.sqrtA = str2num_str(buff[1].substr(61, 19));
                    nav.toe.sec = str2num_str(buff[2].substr(4, 19));
                    nav.toe.week = nav.toc.week;
                    nav.cic = str2num_str(buff[2].substr(23, 19));
                    nav.Omega0 = str2num_str(buff[2].substr(42, 19));
                    nav.cis = str2num_str(buff[2].substr(61, 19));
                    nav.i0 = str2num_str(buff[3].substr(4, 19));
                    nav.crc = str2num_str(buff[3].substr(23, 19));
                    nav.omega = str2num_str(buff[3].substr(42, 19));
                    nav.OmegaDot = str2num_str(buff[3].substr(61, 19));
                    nav.IDot = str2num_str(buff[4].substr(4, 19));
                    if (buff.size() >= 6) nav.tgd = str2num_str(buff[5].substr(42, 19));
                }
            }
        }
    }
    if (nav.sat != 0) navList.push_back(nav);
    std::cout << "已加载 " << navList.size() << " 条导航记录。" << std::endl;
    return navList;
}
// 确定观测文件中的信号列索引
static ColIdx map_signal_col(char sys, const std::vector<std::string>& types) {
    ColIdx idx;
    for (const std::string& target : prio_P1[sys]) {
        auto it = find(types.begin(), types.end(), target);
        if (it != types.end()) { idx.p1 = distance(types.begin(), it); idx.c1 = target; break; }
    }
    if (idx.p1 == -1) {
        for (size_t i = 0; i < types.size(); i++) {
            if (types[i][0] != 'C') continue;
            if ((sys == 'C' && (types[i][1] == '2' || types[i][1] == '1')) || (sys != 'C' && types[i][1] == '1')) {
                idx.p1 = i; idx.c1 = types[i]; break;
            }
        }
    }
    if (idx.p1 != -1) {
        std::string tD = idx.c1; tD[0] = 'D';
        auto it = find(types.begin(), types.end(), tD);
        if (it != types.end()) idx.d1 = distance(types.begin(), it);
    }
    for (const std::string& target : prio_P2[sys]) {
        auto it = find(types.begin(), types.end(), target);
        if (it != types.end()) { idx.p2 = distance(types.begin(), it); idx.c2 = target; break; }
    }
    if (idx.p2 == -1) {
        for (size_t i = 0; i < types.size(); i++) {
            if ((int)i == idx.p1) continue;
            if (types[i][0] != 'C') continue;
            if ((sys == 'G' || sys == 'R') && types[i][1] == '2') { idx.p2 = i; idx.c2 = types[i]; break; }
            if (sys == 'E' && (types[i][1] == '5' || types[i][1] == '7')) { idx.p2 = i; idx.c2 = types[i]; break; }
            if (sys == 'C' && (types[i][1] == '6' || types[i][1] == '7')) { idx.p2 = i; idx.c2 = types[i]; break; }
        }
    }
    if (idx.p2 != -1) {
        std::string tD = idx.c2; tD[0] = 'D';
        auto it = find(types.begin(), types.end(), tD);
        if (it != types.end()) idx.d2 = distance(types.begin(), it);
    }
    return idx;
}

// 读取观测文件
std::vector<std::vector<ObsData>> readObsFile(const std::string& filename, double* approxPos) {
    std::vector<std::vector<ObsData>> allEpochs;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法打开观测文件 " << filename << std::endl;
        return allEpochs;
    }
    std::string line;
    std::map<char, ColIdx> sys_col_idx;
    while (getline(file, line)) {
        if (line.length() < 60) continue;
        if (line.find("APPROX POSITION XYZ") != std::string::npos) {
            approxPos[0] = str2num_str(line.substr(0, 14));
            approxPos[1] = str2num_str(line.substr(14, 14));
            approxPos[2] = str2num_str(line.substr(28, 14));
        }
        if (line.find("SYS / # / OBS TYPES") != std::string::npos) {
            char sys = line[0];
            if (sys == ' ') continue;
            try {
                int n_obs = stoi(line.substr(3, 3));
                std::vector<std::string> types;
                for (int i = 0; i < 13 && i < n_obs; i++) {
                    if (line.length() >= 7 + i * 4 + 3) types.push_back(line.substr(7 + i * 4, 3));
                }
                int lines_to_read = (n_obs - 1) / 13;
                for (int k = 0; k < lines_to_read; k++) {
                    getline(file, line);
                    for (int i = 0; i < 13; i++) {
                        int idx = (k + 1) * 13 + i;
                        if (idx >= n_obs) break;
                        if (line.length() >= 7 + i * 4 + 3) types.push_back(line.substr(7 + i * 4, 3));
                    }
                }
                sys_col_idx[sys] = map_signal_col(sys, types);
             //   std::cout << "[RINEX] " << sys << " 信号映射: P1=" << sys_col_idx[sys].c1 << " P2=" << sys_col_idx[sys].c2 << std::endl;
            }
            catch (...) {}
        }
        if (line.find("END OF HEADER") != std::string::npos) break;
    }
    std::vector<ObsData> epochObs;
    GPSTime epochTime = { 0, 0 };
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!epochObs.empty()) { allEpochs.push_back(epochObs); epochObs.clear(); }
            try {
                int yr = stoi(line.substr(2, 4));
                int mon = stoi(line.substr(7, 2));
                int day = stoi(line.substr(10, 2));
                int hr = stoi(line.substr(13, 2));
                int min = stoi(line.substr(16, 2));
                double sec = str2num_str(line.substr(19, 11));
                double ep[6] = { (double)yr, (double)mon, (double)day, (double)hr, (double)min, sec };
                epochTime = epoch2time(ep);
            }
            catch (...) { continue; }
        }
        else {
            char sys = line[0];
            if (sys_col_idx.find(sys) == sys_col_idx.end()) continue;
            try {
                ObsData obs = { 0 };
                obs.time = epochTime; obs.sys = sys;
                obs.sat = stoi(line.substr(1, 2));
                ColIdx idx = sys_col_idx[sys];
                auto read_val = [&](int col) -> double {
                    if (col < 0) return 0.0;
                    int start = 3 + col * 16;
                    if (line.length() < start + 14) return 0.0;
                    return str2num_str(line.substr(start, 14));
                    };
                obs.P1 = read_val(idx.p1); obs.P2 = read_val(idx.p2);
                // D1 and D2 reading logic kept to maintain file pointer stability, but not used
                // obs.D1 = read_val(idx.d1); obs.D2 = read_val(idx.d2);
                if (obs.P1 != 0.0) epochObs.push_back(obs); // Modified check: only P1 matters now
            }
            catch (...) {}
        }
    }
    if (!epochObs.empty()) allEpochs.push_back(epochObs);
    std::cout << "已加载 " << allEpochs.size() << " 个历元。" << std::endl;
    return allEpochs;
}

// 工具函数 
void ecef2pos_common(const double* r, double* pos) {// ECEF (XYZ) 转 大地坐标
    double e2 = FE_WGS84 * (2.0 - FE_WGS84);
    double r2 = r[0] * r[0] + r[1] * r[1];
    double z = r[2];
    double zk = 0.0;
    double v = RE_WGS84;
    double sinp = 0.0;
    for (double z_old = -1e16; fabs(zk - z_old) > 1e-4;) {
        z_old = zk;
        sinp = z / sqrt(r2 + z * z);
        v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);
        zk = r[2] + v * e2 * sinp;
    }
    pos[0] = (r2 > 1e-12) ? atan(zk / sqrt(r2)) : (r[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
    pos[1] = (r2 > 1e-12) ? atan2(r[1], r[0]) : 0.0;
    pos[2] = sqrt(r2 + z * z) - v;
}

void xyz2enu_common(const double* r, const double* ref, double* enu) {// XYZ 转 站心地平坐标系
    double lat = atan2(ref[2], sqrt(ref[0] * ref[0] + ref[1] * ref[1]));
    double lon = atan2(ref[1], ref[0]);
    double sinP = sin(lat), cosP = cos(lat);
    double sinL = sin(lon), cosL = cos(lon);
    double dx = r[0] - ref[0];
    double dy = r[1] - ref[1];
    double dz = r[2] - ref[2];
    enu[0] = -sinL * dx + cosL * dy;
    enu[1] = -sinP * cosL * dx - sinP * sinL * dy + cosP * dz;
    enu[2] = cosP * cosL * dx + cosP * sinL * dy + sinP * dz;
}

// 计算卫星的方位角和高度角
void sat_azel_common(const double* r, const double* rs, double* azel) {
    double enu[3];
    xyz2enu_common(rs, r, enu);
    double dist = sqrt(enu[0] * enu[0] + enu[1] * enu[1] + enu[2] * enu[2]);
    if (dist < 1e-3) { azel[0] = 0; azel[1] = 0; return; }
    azel[1] = asin(enu[2] / dist);
    azel[0] = atan2(enu[0], enu[1]);
    if (azel[0] < 0) azel[0] += 2 * PI;
}

// 模块实现: 多系统 SPP 核心！
const bool USE_GPS = 1;
const bool USE_BDS = 1;
const bool USE_GAL = 1;
const bool USE_GLO = 0; // 没用GLONASS

static bool is_sys_enabled(char sys) {
    if (sys == 'G') return USE_GPS;
    if (sys == 'C') return USE_BDS;
    if (sys == 'E') return USE_GAL;
    if (sys == 'R') return USE_GLO;
    return false;
}

// 根据系统获取载波频率 用于无电离层组合
static bool get_sys_freqs(char sys, int k, double& f1, double& f2) {
    switch (sys) {
    case 'G': f1 = FREQ_GPS_L1; f2 = FREQ_GPS_L2; return true;
    case 'C': f1 = FREQ_BDS_B1; f2 = FREQ_BDS_B3; return true;
    case 'E': f1 = FREQ_GAL_E1; f2 = FREQ_GAL_E5a; return true;
    case 'R':
        f1 = FREQ_GLO_G1_BASE + k * FREQ_GLO_G1_DELTA;
        f2 = FREQ_GLO_G2_BASE + k * FREQ_GLO_G2_DELTA;
        return true;
    default: return false;
    }
}

// 计算无电离层组合观测值
static double calc_IF(double obs1, double obs2, double f1, double f2) {
    double f1_sq = f1 * f1;
    double f2_sq = f2 * f2;
    return (f1_sq * obs1 - f2_sq * obs2) / (f1_sq - f2_sq);
}

// 单历元定位主解算函数
// 输入: 当前历元观测值, 星历
// 输入输出: 接收机位置, 接收机钟差
// 输出: PDOP值
static int estPosition_GNSS(const std::vector<ObsData>& obsList, const std::vector<NavData>& navList,
    double* rr, std::map<char, double>& clk, double* pdop) {
    std::vector<char> syss;
    for (const auto& o : obsList) {
        if (!is_sys_enabled(o.sys)) continue;
        if (o.P1 != 0 && o.P2 != 0) {
            bool found = false; for (char s : syss) if (s == o.sys) found = true;
            if (!found) syss.push_back(o.sys);
        }
    }
    if (syss.empty()) return 0;
    int m = 3 + syss.size();
    for (char s : syss) if (clk.find(s) == clk.end()) clk[s] = 0.0;
    int nSat = 0;
    for (int iter = 0; iter < 10; iter++) {// 迭代最小二乘 
        std::vector<double> H, v, W;  // H:设计矩阵, v:残差向量, W:权矩阵
        nSat = 0;
        bool valid_pos = (sqrt(rr[0] * rr[0] + rr[1] * rr[1] + rr[2] * rr[2]) > 1.0);
        // 如果当前位置 rr 有效，转为 BLH 以计算对流层
        double pos_blh[3] = { 0 };
        if (valid_pos) ecef2pos_common(rr, pos_blh);
        for (const auto& obs : obsList) {    // 遍历每颗观测卫星
            if (!is_sys_enabled(obs.sys)) continue;
            if (obs.P1 == 0 || obs.P2 == 0) continue;
            const NavData* nav = nullptr; double min_dt = 1e9;
            for (const auto& n : navList) {
                if (n.sat == obs.sat && n.sys == obs.sys) {
                    double dt = timediff(obs.time, n.toe);
                    if (fabs(dt) < fabs(min_dt)) { min_dt = dt; nav = &n; }
                }
            }
            if (!nav || fabs(min_dt) > 14400) continue;
            double f1, f2;
            if (!get_sys_freqs(obs.sys, nav->freq_num, f1, f2)) continue;
            double tau = obs.P1 / CLIGHT;
            GPSTime t_tx = obs.time; t_tx.sec -= tau;
            double rs[3] = { 0 }, dts = 0.0;
            satPos(t_tx, *nav, rs, &dts);
            if (obs.sys == 'C') {
                double gamma = (f1 * f1) / (f1 * f1 - f2 * f2);
                dts -= nav->tgd * gamma;
            }
            double theta = OMEGA_E * tau;
            double rx = rs[0] * cos(theta) + rs[1] * sin(theta);
            double ry = -rs[0] * sin(theta) + rs[1] * cos(theta);
            rs[0] = rx; rs[1] = ry;
            double r = sqrt(pow(rs[0] - rr[0], 2) + pow(rs[1] - rr[1], 2) + pow(rs[2] - rr[2], 2));
            double azel[2] = { 0 };
            if (valid_pos) {
                sat_azel_common(rr, rs, azel);
                if (azel[1] < EL_MASK) continue;
            }
            else {
                azel[1] = PI / 2.0;
            }
            double trop = 0, var_trop = 0;
            if (iter > 0 && valid_pos) {
                trop_model_prec(obs.time, pos_blh, azel, &trop, &var_trop);
            }
            double P_IF = calc_IF(obs.P1, obs.P2, f1, f2);
            double res = P_IF - (r + clk[obs.sys] - CLIGHT * dts + trop);
            H.push_back((rr[0] - rs[0]) / r);
            H.push_back((rr[1] - rs[1]) / r);
            H.push_back((rr[2] - rs[2]) / r);
            for (char s : syss) H.push_back(s == obs.sys ? 1.0 : 0.0);
            v.push_back(res);
            double w = 1.0;
            if (valid_pos) {
                double sin_el = sin(azel[1]);
                w = sin_el * sin_el;
                if (obs.sys == 'C' && obs.sat <= 5) w *= 0.1;
            }
            W.push_back(w);
            nSat++;
        }
        if (nSat < m) return 0;
        std::vector<double> Q(m * m, 0), b(m, 0), Q_geo(m * m, 0);
        for (int i = 0; i < nSat; i++) {
            for (int j = 0; j < m; j++) {
                for (int k = 0; k < m; k++) {
                    Q[j * m + k] += H[i * m + j] * W[i] * H[i * m + k];
                    Q_geo[j * m + k] += H[i * m + j] * H[i * m + k];
                }
                b[j] += H[i * m + j] * W[i] * v[i];
            }
        }
        if (matinv(Q.data(), m) == -1) return 0;
        std::vector<double> dx(m);
        matmul("NN", m, m, 1, 1.0, Q.data(), b.data(), 0.0, dx.data());
        rr[0] += dx[0]; rr[1] += dx[1]; rr[2] += dx[2];
        for (int k = 0; k < (int)syss.size(); k++) clk[syss[k]] += dx[3 + k];
        if (matinv(Q_geo.data(), m) != -1) *pdop = sqrt(Q_geo[0] + Q_geo[m + 1] + Q_geo[2 * m + 2]);
        if (sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) < 1e-4) return nSat;
    }
    return nSat;
}

// 进行批处理
void spp_gnss_process(const std::vector<std::vector<ObsData>>& obsList,
    const std::vector<NavData>& navList,
    const double* refPos,
    const std::string& outfile)
{
    
    ofstream fout(outfile);
    if (!fout.is_open()) {
        cerr << "无法打开输出文件: " << outfile << endl;
        return;
    }

    fout << "多系统双频无电离层组合单点定位 " << endl;
    fout << "启用系统: ";
    if (USE_GPS) fout << "GPS ";
    if (USE_BDS) fout << "BDS ";
    if (USE_GAL) fout << "GAL ";
    if (USE_GLO) fout << "GLO ";
    fout << endl;

    fout << fixed << setprecision(3);
    fout << "Time(GPST)  卫星数  PDOP      dE(米)    dN(米)   dU(米)" << endl;

    for (const auto& epochObs : obsList) {
        if (epochObs.empty()) continue;
        GPSTime t = epochObs[0].time;
        double ep[6]; time2epoch(t, ep);
        double current_hour = ep[3] + ep[4] / 60.0 + ep[5] / 3600.0;

        //  02:00 - 06:00
        if (current_hour < 2.0 - 1e-9 || current_hour > 6.0 + 1e-9) continue;

        double rr[3] = { 0 };
        if (refPos[0] != 0) { rr[0] = refPos[0]; rr[1] = refPos[1]; rr[2] = refPos[2]; }
        std::map<char, double> clk;
        double pdop = 0.0;
        int ns = estPosition_GNSS(epochObs, navList, rr, clk, &pdop);
        if (ns >= 4) {

            if (ep[5] >= 59.95) {
                ep[5] = 0; ep[4]++;
                if (ep[4] >= 60) { ep[4] = 0; ep[3]++; }
            }
            char time_str[32];
            sprintf(time_str, "%02.0f:%02.0f:%04.1f", ep[3], ep[4], ep[5]);
            double enu_pos[3] = { 0 };
            if (refPos[0] != 0) xyz2enu_common(rr, refPos, enu_pos);
            else xyz2enu_common(rr, rr, enu_pos);

            fout << time_str << "   "
                << setw(2) << ns << "   "
                << setw(5) << pdop << "   "
                << setw(8) << enu_pos[0] << " "
                << setw(8) << enu_pos[1] << " "
                << setw(8) << enu_pos[2] << endl;
        }
    }
    cout << "处理完成。结果已写入 " << outfile << endl;
    fout.close();
}


int main() {
    string nav_file = "brdc1590.24p";
    string obs_file = "jfng1590.24o";
    string out_file = "spp_result.txt";

    cerr << "GNSS (GPS/GAL/BDS) 标准单点定位程序 " << endl;
    cerr << "正在读取导航电文文件: " << nav_file << endl;
    auto navList = readNavFile(nav_file);
    if (navList.empty()) {
        cerr << "读取导航文件失败。" << endl;
        return 1;
    }

    cerr << "正在读取观测文件: " << obs_file << endl;
    double refPos[3] = { 0 };
    auto obsList = readObsFile(obs_file, refPos);
    if (obsList.empty()) {
        cerr << "读取观测文件失败。" << endl;
        return 1;
    }

    spp_gnss_process(obsList, navList, refPos, out_file);

    return 0;
}