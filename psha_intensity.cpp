#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <random>
#include <omp.h>
#include <iomanip>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ================= КОНСТАНТЫ ИЗ ОРИГИНАЛА =================
static const double RAD = 180.0 / M_PI;  // Коэффициент перевода радиан в градусы
static const double GI0_DEFAULT = 4.5;   // Начальная интенсивность гистограммы
static const double DI = 0.1;            // Шаг гистограммы
static const int IMGS = 80;              // Размер гистограммы
static const double R3DMIN = 3.0;        // Минимальное расстояние

// ================= МОДЕЛЬ ЗАТУХАНИЯ =================
struct AttenModel {
    std::string name;
    double dida;      // Коэффициент при log10(att)
    double didmw;     // Коэффициент при (Mw - Mw_bas)
    double cmag;      // Базовая магнитуда для площади разрыва
    double an1, rq1;  // Параметры затухания для ближней зоны
    double rswitch;   // Точка переключения
    double an2, rq2;  // Параметры затухания для дальней зоны
    double amwbas;    // Базовая магнитуда
    double rbas;      // Базовое расстояние
    double aibas;     // Базовая интенсивность
    double albywb;    // L/W ratio для базовой магнитуды
    double gi0;       // Начало гистограммы
    // CLW1, CLW2, CMW1, CMW2 для интерполяции L/W
    double clw1, clw2, cmw1, cmw2;
    // r35..r39 для радиусов влияния
    double r35, r36, r37, r38, r39;
};

// Встроенные модели (значения из оригинала)
std::vector<AttenModel> builtin_models = {
    {"45imr", 1.667, 1.85, 4.1, 1.0, 100, 70, 0.5, 100, 6.23, 50, 6.05, 2.0, 4.5, 2.0, 2.0, 4.0, 8.0, 105, 270, 470, 820, 1000},
    {"11imr", 1.667, 2.15, 4.1, 1.2, 150, 500, 1.1, 150, 6.97, 117, 6.00, 3.0, 4.5, 2.0, 3.0, 4.0, 8.0, 70, 150, 240, 320, 400},
    {"41imr", 1.667, 1.85, 4.1, 1.0, 90, 1000, 1.0, 90, 8.00, 100, 7.70, 4.0, 4.5, 2.0, 4.0, 4.0, 8.0, 105, 270, 470, 820, 1000},
    {"kim2",  1.667, 1.85, 4.1, 0.9, 100, 5000, 0.5, 100, 6.97, 130, 6.00, 2.0, 4.5, 2.0, 2.0, 4.0, 8.0, 105, 270, 470, 820, 1000},
    {"KLN",   1.667, 1.85, 4.1, 0.58, 100, 70, 0.5, 100, 4.80, 55, 4.50, 2.0, 4.5, 2.0, 2.0, 4.0, 8.0, 105, 270, 470, 820, 1000},
};

const AttenModel* findModel(const std::string& name) {
    for (auto& m : builtin_models)
        if (m.name == name) return &m;
    return nullptr;
}

// ================= СОБЫТИЕ ИЗ CTL =================
struct Event {
    int ind;          // Индекс зоны
    double mag;       // Магнитуда Mw
    double L, W;      // Длина и ширина разрыва (км)
    double az;        // Азимут (градусы)
    double dip;       // Угол падения (градусы)
    double lat, lon;  // Координаты угла X1 (градусы)
    double depth;     // Глубина (км)
};

// ================= ПАРАМЕТРЫ ЗОНЫ =================
struct ZoneParams {
    //double akdip, akaz;       // Смещение центра относительно площадки
    double sdevm, sdevi;      // СКО для магнитудной и индивидуальной поправки
    std::string model_name;   // Имя модели затухания
    //int tip_podv;             // Тип подвижки (не используется в расчете интенсивности)
    double devl, devc;        // СКО для смещения эпицентра
    double crot, srot;        // cos/sin азимута (вычисляются)
};

// ================= ГЕОДЕЗИЧЕСКИЕ ПРЕОБРАЗОВАНИЯ (ТОЧНАЯ КОПИЯ ИЗ C#) =================

// GEDECCON: Преобразование географических координат в плоские (км)
// f0, al0, az0 - координаты и азимут центра проекции (радианы)
// f, al - географические координаты точки (радианы)
// x, y - выходные плоские координаты (км)
void GEDECCON(double f0, double al0, double az0, double f, double al, double& x, double& y) {
    const double num1 = 10374.71;
    const double num2 = 0.854116;
    
    double num3 = num1 - (6367.5584958746 * f0 - 16.0364802690885 * std::sin(2.0 * f0) 
                  + 0.0168280667831 * std::sin(4.0 * f0) - 2.1975279E-05 * std::sin(6.0 * f0) 
                  + 3.11243E-08 * std::sin(8.0 * f0));
    
    double num4 = num1 - (6367.5584958746 * f - 16.0364802690885 * std::sin(2.0 * f) 
                  + 0.0168280667831 * std::sin(4.0 * f) - 2.1975279E-05 * std::sin(6.0 * f) 
                  + 3.11243E-08 * std::sin(8.0 * f));
    
    double num5 = num2 * (al - al0);
    double num6 = num4 * std::sin(num5);
    double num7 = num3 - num4 * std::cos(num5);
    double num8 = std::sin(az0);
    double num9 = std::cos(az0);
    
    x = num6 * num9 - num7 * num8;
    y = num6 * num8 + num7 * num9;
}

// DEGEDCON: Обратное преобразование плоских координат в географические
void DEGEDCON(double f0, double al0, double az0, double x, double y, double& f, double& al) {
    const double num1 = 6367.5584958746;
    const double num2 = 10374.71;
    const double num3 = 0.854116;
    
    double num4 = num2 - (6367.5584958746 * f0 - 16.0364802690885 * std::sin(2.0 * f0) 
                  + 0.0168280667831 * std::sin(4.0 * f0) - 2.1975279E-05 * std::sin(6.0 * f0) 
                  + 3.11243E-08 * std::sin(8.0 * f0));
    
    double num5 = x * std::cos(az0) + y * std::sin(az0);
    double num6 = -x * std::sin(az0) + y * std::cos(az0);
    double num7 = std::sqrt((num4 - num6) * (num4 - num6) + num5 * num5);
    
    f = (num2 - num7) / num1 + 0.0024412912;
    double num8 = std::asin(num5 / num7);
    al = num8 / num3 + al0;
}

// ================= ФУНКЦИИ ЗАТУХАНИЯ (ТОЧНАЯ КОПИЯ ИЗ C#) =================

// CATT: Базовая функция затухания
inline double CATT(double r, double an, double rq) {
    double r_pow_an = std::pow(r, an);
    return 1.0 / (r_pow_an * r_pow_an * std::exp(r / rq));
}

// ATT: Полная функция затухания с переключением
double ATT(double r, double /*mag*/, const AttenModel& m) {
    double RA = m.rswitch * 0.85;
    double RB = m.rswitch * 1.15;
    double DR = RB - RA;
    
    // Вычисляем CC (коэффициент сшивки)
    double CC = CATT(m.rswitch, m.an1, m.rq1) / CATT(m.rswitch, m.an2, m.rq2);
    
    if (r <= RA) {
        return CATT(r, m.an1, m.rq1);
    } else if (r >= RB) {
        return CC * CATT(r, m.an2, m.rq2);
    } else {
        // Линейная интерполяция в переходной зоне
        double att_ra = CATT(RA, m.an1, m.rq1);
        double att_rb = CC * CATT(RB, m.an2, m.rq2);
        return att_ra + (r - RA) / DR * (att_rb - att_ra);
    }
}

// ================= FINCOR: ПОПРАВКА ЗА КОНЕЧНОСТЬ ИСТОЧНИКА (ТОЧНАЯ КОПИЯ) =================

// Вычисляет поправку FCOR и минимальное расстояние AMINDIST
// r - горизонтальное расстояние до приёмника
// SBEt, CBEt - sin/cos угла между направлением на приёмник и осью X
// RPAR - параметры площадки [1]=L, [2]=W, [3]=H (глубина центра), [4]=phi (угол наклона, =DIP-PI/2), [5]=NL, [6]=NW
void FINCOR(double r, double SBEt, double CBEt, double& AMINDIST, double& corr, 
            double mag, const AttenModel& m, const double* RPAR) {
    double AL = RPAR[1];   // Длина (по простиранию)
    double AW = RPAR[2];   // Ширина (по падению)
    double H = RPAR[3];    // Глубина центра (положительное значение = вниз)
    double phi = RPAR[4];  // Угол = DIP - PI/2 (отрицательный для нормального падения)
    double NL = RPAR[5];   // Число точек по длине
    double NW = RPAR[6];   // Число точек по ширине
    
    double NL2 = NL / 2.0;
    double NW2 = NW / 2.0;
    double DL = AL / NL;
    double DW = AW / NW;
    
    // Координаты приёмника (горизонтальные, относительно центра площадки)
    double XRC = r * CBEt;  // По X (перпендикулярно простиранию)
    double YRC = r * SBEt;  // По Y (вдоль простирания)
    
    // Cos и sin угла phi (= DIP - PI/2)
    // Если DIP = 90° (вертикальный), то phi = 0, CP = 1, SP = 0
    // Если DIP = 0° (горизонтальный), то phi = -PI/2, CP = 0, SP = -1
    double CP = std::cos(phi);
    double SP = std::sin(phi);
    
    AMINDIST = 1e7;
    double SUM = 0.0;
    double SUMW = 0.0;
    
    // Обход точек на площадке разрыва
    for (int i = 1; i <= (int)NL; ++i) {
        for (int j = 1; j <= (int)NW; ++j) {
            // YS - координата по простиранию (горизонтальная, "вдоль" разлома)
            double YS = DL * ((double)(i - 1) - NL2);
            
            // XS, ZS - координаты по падению
            // XS - горизонтальная проекция (перпендикулярно простиранию)
            // ZS - вертикальная компонента (отрицательная = вниз)
            double XS = DW * ((double)(j - 1) - NW2) * SP;
            double ZS = DW * ((double)(j - 1) - NW2) * CP - H;  // -H потому что глубина положительная
            
            // Расстояние от точки на площадке до приёмника
            double DIST = std::sqrt(std::pow(YS - YRC, 2.0) + std::pow(XS - XRC, 2.0) + std::pow(ZS, 2.0));
            
            double CONTRIB = ATT(DIST, mag, m);
            SUM += CONTRIB;
            SUMW += 1.0;
            
            if (DIST < AMINDIST) AMINDIST = DIST;
        }
    }
    
    // Расстояние от приёмника до центра (3D)
    double RR = std::sqrt(r * r + H * H);
    double att_center = ATT(RR, mag, m);
    
    corr = (SUM / SUMW) / att_center;
}

// ================= ВЫЧИСЛЕНИЕ ALBYW (L/W RATIO) =================

double computeALBYW(double mag, const AttenModel& m) {
    if (mag < m.cmw1) return m.clw1;
    if (mag >= m.cmw2) return m.clw2;
    return m.clw1 + (m.clw2 - m.clw1) * (mag - m.cmw1) / (m.cmw2 - m.cmw1);
}

// ================= ПРЕДВЫЧИСЛЕНИЕ DENOM =================

double computeDenom(const AttenModel& m) {
    // Из MACRR3 (LLOOP == 1):
    // SB = 10^(AMWBAS - CMAG)
    // ALBYWB = CLW для AMWBAS (интерполяция)
    // ALB = sqrt(ALBYWB * SB)
    // AWB = SB / ALB
    
    double SB = std::pow(10.0, m.amwbas - m.cmag);
    double ALBYWB = computeALBYW(m.amwbas, m);
    double ALB = std::sqrt(ALBYWB * SB);
    double AWB = SB / ALB;
    
    // Дискретизация (шаг 2 км)
    double DL = 2.0;
    double DW = DL;
    
    double NLB = ALB / DL;
    if (NLB > 199.0) NLB = 301.0;
    NLB = (int)(NLB / 2.0) * 2.0 + 1.0;
    
    double NWB = AWB / DW;
    if (NWB > 99.0) NWB = 201.0;
    NWB = (int)(NWB / 2.0) * 2.0 + 1.0;
    
    double HB = 0.0;   // Глубина центра базового источника
    double PHIB = 0.0; // Угол падения базового источника (DIP - PI/2, для DIP=90 -> PHI=0)
    
    // Заполняем RPAR
    double RPAR[7];
    RPAR[1] = ALB;
    RPAR[2] = AWB;
    RPAR[3] = HB;
    RPAR[4] = PHIB;
    RPAR[5] = NLB;
    RPAR[6] = NWB;
    
    double DISTMIN, FCORB;
    // В оригинале: FINCOR(RBAS, SBEt=0, CBEt=1, ...)
    // Приёмник на расстоянии RBAS по оси X (перпендикулярно простиранию)
    FINCOR(m.rbas, 0.0, 1.0, DISTMIN, FCORB, m.amwbas, m, RPAR);
    
    // DENOM = FCORB * ATT(RBAS, AMWBAS)
    return FCORB * ATT(m.rbas, m.amwbas, m);
}

// ================= ВЫЧИСЛЕНИЕ РАДИУСА ВЛИЯНИЯ =================

double computeRBALL3(double mag, const AttenModel& m) {
    // Таблица BSM3, DST3 из оригинала
    static const double BSM3[] = {0, 0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double DST3[] = {0, 0, m.r35/2.0, m.r35, m.r36, m.r37, m.r38, m.r39};
    
    const int IM3 = 7;
    double r3m = DST3[IM3];
    
    for (int i = 1; i < IM3; ++i) {
        if (mag >= BSM3[i] && mag <= BSM3[i+1]) {
            double t = (mag - BSM3[i]) / (BSM3[i+1] - BSM3[i]);
            r3m = DST3[i] + t * (DST3[i+1] - DST3[i]);
            break;
        }
    }
    
    return r3m;
}

// ================= СТРУКТУРА ДЛЯ КЭШИРОВАННОГО СОБЫТИЯ =================

struct CachedEvent {
    // Плоские координаты центра (XC) относительно центра проекции
    double XC_x, XC_y, XC_z;  // SPAR[9], SPAR[10], SPAR[11]
    
    // Параметры источника
    double mag;       // SPAR[1]
    double SL0, SW0;  // SPAR[2], SPAR[3] - размеры
    double SRCAZ;     // SPAR[4] - азимут (радианы, в локальной системе)
    double DIP;       // SPAR[5] - угол падения (радианы)
    
    // Косинус и синус азимута
    double CROT, SROT;
    
    // Параметры для RPAR
    double RPAR[7];
    
    // Модель затухания
    const AttenModel* model;
    double DENOM;
    
    // Случайные поправки (фиксированные для данного события)
    double UI, VI;      // Для индивидуальной поправки
    double GLBVI;       // = SDEVM * VI
    
    // Параметры зоны
    double SDEVI, SDEVM;
    
    // Радиус влияния
    double RBALL3;
};

// ================= ЗАГРУЗКА ПАРАМЕТРОВ ЗОН =================

std::map<int, ZoneParams> loadZoneParams(const std::string& fname) {
    std::map<int, ZoneParams> params;
    std::ifstream fin(fname);
    if (!fin) {
        std::cerr << "Error: cannot open " << fname << std::endl;
        exit(1);
    }
    std::string line;
    std::getline(fin, line); // header
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        int ind;
        ZoneParams zp;
        if (!(iss >> ind  >> zp.sdevm >> zp.sdevi 
              >> zp.model_name >> zp.devl >> zp.devc)) {
            // Если devl, devc не указаны, используем 0
            zp.devl = 0.0;
            zp.devc = 0.0;
        }
        params[ind] = zp;
    }
    return params;
}

// ================= ГЛАВНЫЙ КЛАСС РАСЧЕТА =================

class IntensityCalculator {
public:
    // Центр проекции
    double PHI0, AL0, AZ0;
    
    // Кэш DENOM для каждой модели
    std::map<std::string, double> denom_cache;
    
    // Генератор случайных чисел
    std::mt19937 rng;
    std::normal_distribution<double> norm_dist;
    
    IntensityCalculator() : rng(345), norm_dist(0.0, 1.0) {
        PHI0 = AL0 = AZ0 = 0.0;
        
        // Предвычисляем DENOM для всех моделей
        for (const auto& m : builtin_models) {
            denom_cache[m.name] = computeDenom(m);
        }
    }
    
    // Установка центра проекции по набору точек
    void setProjectionCenter(const std::vector<std::pair<double,double>>& points) {
        if (points.empty()) return;
        
        double sum_lat = 0, sum_lon = 0;
        for (const auto& p : points) {
            sum_lat += p.first;
            sum_lon += p.second;
        }
        
        PHI0 = (sum_lat / points.size()) / RAD;  // в радианы
        AL0 = (sum_lon / points.size()) / RAD;
        AZ0 = 0.0;
    }
    
    // Преобразование события из CTL в кэшированное
    CachedEvent cacheEvent(const Event& e, const ZoneParams& zp) {
        CachedEvent ce;
        
        // Находим модель
        ce.model = findModel(zp.model_name);
        if (!ce.model) {
            std::cerr << "Unknown model: " << zp.model_name << std::endl;
            exit(1);
        }
        ce.DENOM = denom_cache[zp.model_name];
        
        // Магнитуда и размеры
        ce.mag = e.mag;
        ce.SL0 = e.L;
        ce.SW0 = e.W;
        
        // Преобразуем координаты угла X1 в плоские
        double X1_lat_rad = e.lat / RAD;
        double X1_lon_rad = e.lon / RAD;
        double X1_x, X1_y;
        GEDECCON(PHI0, AL0, AZ0, X1_lat_rad, X1_lon_rad, X1_x, X1_y);
        double X1_z = e.depth;  // Глубина
        
        // Азимут и угол падения (переводим в радианы и в локальную систему)
        ce.SRCAZ = e.az / RAD - AZ0;
        ce.DIP = e.dip / RAD;
        
        ce.CROT = std::cos(ce.SRCAZ);
        ce.SROT = std::sin(ce.SRCAZ);
        
        // Вычисляем центр XC из угла X1
        // В оригинале: SRCR вычисляет X1,X2,X3,X4,XC из X0 (гипоцентра)
        // Но в CTL записывается X1 (угол), поэтому нужно обратное преобразование
        
        double cos_az = std::cos(ce.SRCAZ);
        double sin_az = std::sin(ce.SRCAZ);
        double cos_dip = std::cos(ce.DIP);
        double sin_dip = std::sin(ce.DIP);
        
        // Единичные векторы по простиранию и падению
        double eL_x = sin_az;
        double eL_y = cos_az;
        double eL_z = 0.0;
        
        double eW_x = cos_az * cos_dip;
        double eW_y = -sin_az * cos_dip;
        double eW_z = sin_dip;
        
        // Центр относительно X1 (угол -L/2-akaz*L по простиранию, -W/2-akdip*W по падению)
        // На самом деле в SRCR:
        // X1 = XC - 0.5*W*eW - 0.5*L*eL + L*eL (последнее - сдвиг, см. код)
        // Упрощённо: XC = X1 + смещения
        
        // Из кода SRCR (C#):
        // num4 = -AKDIP * SW0
        // num5 = -AKAZ * SL0
        // XC[index] = X0[index] + num4 * eW[index] + num5 * eL[index]
        // num6 = 0.5 * SL0 * eL
        // num7 = 0.5 * SW0 * eW
        // X1[index] = XC[index] - num7 - num6 + num6 = XC[index] - num7
        // 
        // То есть X1 = XC - 0.5*W*eW
        // Следовательно: XC = X1 + 0.5*W*eW
        //
        // НО! В CTL записываются SPAR[6,7,8] = X1, а не X0!
        // А центр SPAR[9,10,11] = XC
        // При этом в основном цикле используется именно XC (SPAR[9], SPAR[10], SPAR[11])
        
        ce.XC_x = X1_x + 0.5 * ce.SW0 * eW_x;
        ce.XC_y = X1_y + 0.5 * ce.SW0 * eW_y;
        ce.XC_z = X1_z + 0.5 * ce.SW0 * eW_z;
        
        // Заполняем RPAR для FINCOR
        ce.RPAR[1] = ce.SL0;
        ce.RPAR[2] = ce.SW0;
        ce.RPAR[3] = ce.XC_z;  // Глубина центра
        ce.RPAR[4] = ce.DIP - M_PI / 2.0;  // Угол для FINCOR
        
        // Дискретизация площадки (из оригинала: шаг 2 км, но ограничения на число точек)
        // В основном цикле: 
        // RPAR[5] = 2.0 * ((int)(CRN * SPAR[2] / 2.0) / 2.0) + 1.0
        // RPAR[6] = 2.0 * (int)(SPAR[3] / 2.0 / 2.0) + 1.0
        // Где CRN = 1 или 2 в зависимости от расстояния
        // Для предвычисления используем CRN = 2 (максимальная детализация)
        
        double CRN = 2.0;
        double NL = 2.0 * (double)((int)(CRN * ce.SL0 / 2.0) / 2) + 1.0;
        double NW = 2.0 * (double)((int)(ce.SW0 / 2.0) / 2) + 1.0;
        
        // Ограничения из оригинала
        if (NL < 1.0) NL = 1.0;
        if (NW < 1.0) NW = 1.0;
        
        ce.RPAR[5] = NL;
        ce.RPAR[6] = NW;
        
        // Параметры зоны
        ce.SDEVI = zp.sdevi;
        ce.SDEVM = zp.sdevm;
        
        // Генерируем случайные поправки
        ce.UI = norm_dist(rng);
        ce.VI = norm_dist(rng);
        ce.GLBVI = ce.SDEVM * ce.VI;
        
        // Радиус влияния
        ce.RBALL3 = computeRBALL3(ce.mag, *ce.model);
        
        return ce;
    }
    
    // Вычисление интенсивности в точке от одного события
    double computeIntensity(const CachedEvent& ce, double site_x, double site_y) {
        // Вектор от центра площадки к точке наблюдения
        double XX = site_x - ce.XC_x;
        double YY = site_y - ce.XC_y;
        
        double RR_sq = XX * XX + YY * YY;
        double RR = std::sqrt(RR_sq);
        
        // 3D расстояние до центра
        double R3D = std::sqrt(RR_sq + ce.XC_z * ce.XC_z);
        
        if (R3D < R3DMIN) R3D = R3DMIN;
        
        // Проверка радиуса влияния
        if (R3D > ce.RBALL3) return -999.0;
        
        if (RR < 1e-5) RR = 0.01;
        
        // В оригинале MACRR3:
        // CBET = (XX * CROT - YY * SROT) / RR
        // SBET = (XX * SROT + YY * CROT) / RR
        // 
        // Это поворот вектора (XX, YY) на угол -SRCAZ
        // CBET = cos(angle_to_receiver - SRCAZ) = cos(angle_in_fault_system)
        // SBET = sin(angle_to_receiver - SRCAZ) = sin(angle_in_fault_system)
        //
        // В системе координат площадки:
        // - ось X направлена по падению (перпендикулярно простиранию)
        // - ось Y направлена по простиранию
        
        double CBET = (XX * ce.CROT - YY * ce.SROT) / RR;
        double SBET = (XX * ce.SROT + YY * ce.CROT) / RR;
        
        // Вычисляем FCOR
        double DISTMIN, FCOR;
        FINCOR(RR, SBET, CBET, DISTMIN, FCOR, ce.mag, *ce.model, ce.RPAR);
        
        // Вычисляем интенсивность по формуле из MACRR3:
        // BALL = AIBAS + DIDMW * (AMW - AMWBAS) + DIDA * log10(FCOR * ATT(R3D) / DENOM)
        double att = ATT(R3D, ce.mag, *ce.model);
        double BALL = ce.model->aibas 
                    + ce.model->didmw * (ce.mag - ce.model->amwbas)
                    + ce.model->dida * std::log10(FCOR * att / ce.DENOM);
        
        // Добавляем случайную поправку:
        // RESI = BALL + UI * SDEVI + GLBVI
        // где GLBVI = SDEVM * VI (глобальная поправка за магнитуду)
        double RESI = BALL + ce.UI * ce.SDEVI + ce.GLBVI;
        
        return RESI;
    }
};

// ================= RSKVN: ВЫЧИСЛЕНИЕ КВАНТИЛЯ ИЗ ГИСТОГРАММЫ =================

// Из оригинала:
// S1[i] = S[i] / PNUM - нормированная кумулятивная гистограмма
// Ищем первый индекс где S1[i] <= 1
// X1 = GI0 + 0.1 * (index1 - 2) + 0.1 * (S1[index1-1] - 1.0) / (S1[index1-1] - S1[index1])
//
// PNUM = NCYCL * (TPR / target_T)
// где TPR - время одного цикла, target_T - период повторяемости (500, 1000, 5000 лет)

double RSKVN(double PNUM, const double* S, double GI0) {
    // S - кумулятивная гистограмма (S[1] - число событий >= GI0+0.0, S[2] - >= GI0+0.1, и т.д.)
    // На самом деле S[1] соответствует GI0, S[2] - GI0+0.1, ...
    
    double SUM = 0.0;
    for (int i = 1; i <= IMGS; ++i) {
        SUM += S[i] * S[i];
    }
    
    if (SUM < 1e-5) {
        return 1.0;  // Нет данных
    }
    
    int index1 = 1;
    for (int i = 1; i <= IMGS; ++i) {
        double S1_i = S[i] / PNUM;
        if (S1_i <= 1.0) {
            index1 = i;
            break;
        }
    }
    
    if (index1 == 1) {
        return 2.0;  // Все значения <= 1, возвращаем минимальную интенсивность
    }
    
    double S1_prev = S[index1 - 1] / PNUM;
    double S1_curr = S[index1] / PNUM;
    
    // Линейная интерполяция
    // X1 = GI0 + 0.1 * (index1 - 2) + 0.1 * (S1_prev - 1.0) / (S1_prev - S1_curr)
    double X1 = GI0 + 0.1 * (double)(index1 - 2) + 0.1 * (S1_prev - 1.0) / (S1_prev - S1_curr);
    
    return X1;
}

// ================= MAIN =================

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " CTL.txt param_table.txt TMAX NCYCL grid_file [output_file]\n";
        std::cerr << "  CTL.txt       - каталог событий\n";
        std::cerr << "  param_table.txt - параметры зон (ind akdip akaz sdevm sdevi parfln tip_podv [devl devc])\n";
        std::cerr << "  TMAX          - время моделирования (лет) для одного цикла\n";
        std::cerr << "  NCYCL         - число циклов\n";
        std::cerr << "  grid_file     - файл с сеткой точек (lat lon)\n";
        std::cerr << "  output_file   - выходной файл (опционально)\n";
        return 1;
    }
    
    std::string ctl_file = argv[1];
    std::string param_file = argv[2];
    double TMAX = std::stod(argv[3]);
    int NCYCL = std::stoi(argv[4]);
    std::string grid_file = argv[5];
    std::string out_file = (argc >= 7) ? argv[6] : "";
    
    double total_years = TMAX * NCYCL;
    
    // Загружаем параметры зон
    auto zoneParams = loadZoneParams(param_file);
    std::cerr << "Loaded " << zoneParams.size() << " zone parameters" << std::endl;
    
    // Читаем события из CTL
    std::vector<Event> events;
    std::ifstream fctl(ctl_file);
    if (!fctl) {
        std::cerr << "Error: cannot open " << ctl_file << std::endl;
        return 1;
    }
    std::string line;
    std::getline(fctl, line); // header
    while (std::getline(fctl, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Event e;
        if (!(iss >> e.ind >> e.mag >> e.L >> e.W >> e.az >> e.dip >> e.lat >> e.lon >> e.depth)) {
            std::cerr << "Warning: bad line in CTL: " << line << std::endl;
            continue;
        }
        events.push_back(e);
    }
    fctl.close();
    std::cerr << "Loaded " << events.size() << " events from CTL" << std::endl;
    
    // Загружаем сетку точек
    std::vector<std::pair<double,double>> grid;
    std::ifstream fgrid(grid_file);
    if (!fgrid) {
        std::cerr << "Error: cannot open grid file " << grid_file << std::endl;
        return 1;
    }
    while (std::getline(fgrid, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        double lat, lon;
        if (!(iss >> lat >> lon)) continue;
        if (std::fabs(lat - 9999.0) < 1e-6) continue;  // Разделитель
        if (std::fabs(lat + 9999.0) < 1e-6) break;     // Конец файла
        grid.push_back({lat, lon});
    }
    fgrid.close();
    std::cerr << "Loaded " << grid.size() << " grid points" << std::endl;
    
    // Создаём калькулятор
    IntensityCalculator calc;
    
    // Устанавливаем центр проекции по сетке
    calc.setProjectionCenter(grid);
    std::cerr << "Projection center: PHI0=" << calc.PHI0 * RAD << ", AL0=" << calc.AL0 * RAD << std::endl;
    
    // Кэшируем события
    std::vector<CachedEvent> cachedEvents;
    cachedEvents.reserve(events.size());
    for (const auto& e : events) {
        auto it = zoneParams.find(e.ind);
        if (it == zoneParams.end()) {
            std::cerr << "Warning: no zone params for ind=" << e.ind << ", using defaults" << std::endl;
            ZoneParams zp;
            zp.sdevm = 0.5;
            zp.sdevi = 0.5;
            zp.model_name = "45imr";
            zp.devl = 0.0;
            zp.devc = 0.0;
            cachedEvents.push_back(calc.cacheEvent(e, zp));
        } else {
            cachedEvents.push_back(calc.cacheEvent(e, it->second));
        }
    }
    std::cerr << "Cached " << cachedEvents.size() << " events" << std::endl;
    
    // Преобразуем координаты сетки в плоские
    std::vector<std::pair<double,double>> grid_xy;
    grid_xy.reserve(grid.size());
    for (const auto& p : grid) {
        double lat_rad = p.first / RAD;
        double lon_rad = p.second / RAD;
        double x, y;
        GEDECCON(calc.PHI0, calc.AL0, calc.AZ0, lat_rad, lon_rad, x, y);
        grid_xy.push_back({x, y});
    }
    
    // Выходной поток
    std::ostream* out = &std::cout;
    std::ofstream fout;
    if (!out_file.empty()) {
        fout.open(out_file);
        if (!fout) {
            std::cerr << "Error: cannot create output file " << out_file << std::endl;
            return 1;
        }
        out = &fout;
    }
    
    // Вычисление T500 для каждой точки
    // PNUM = NCYCL для T = TPR (базовый период)
    // Для T500: PNUM1 = NCYCL * (TPR / 500)
    double target_T = 500.0;  // Период повторяемости в годах
    double PNUM_T500 = (double)NCYCL * (TMAX / target_T);
    
    std::cerr << "TMAX=" << TMAX << ", NCYCL=" << NCYCL << ", PNUM_T500=" << PNUM_T500 << std::endl;
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < grid.size(); ++i) {
        double site_x = grid_xy[i].first;
        double site_y = grid_xy[i].second;
        
        // Накапливаем гистограмму (как в оригинале: IGST[index8, igist])
        std::vector<double> IGST(IMGS + 1, 0.0);
        
        for (const auto& ce : cachedEvents) {
            double resi = calc.computeIntensity(ce, site_x, site_y);
            if (resi < -900) continue;  // Вне радиуса влияния
            
            // Индекс в гистограмме (из оригинала):
            // IGIST = RESI >= GI0 ? (int)((RESI - GI0) / 0.1) + 2 : 1
            int igist;
            if (resi < GI0_DEFAULT) {
                igist = 1;
            } else {
                igist = (int)((resi - GI0_DEFAULT) / DI) + 2;
            }
            if (igist > IMGS) igist = IMGS;
            if (igist < 1) igist = 1;
            
            IGST[igist] += 1.0;
        }
        
        // Преобразуем в кумулятивную гистограмму (из оригинала: S1[i] = S1[i+1] + IGST[i+1])
        std::vector<double> S1(IMGS + 1, 0.0);
        S1[IMGS] = 0.0;
        for (int j = IMGS - 1; j >= 1; --j) {
            S1[j] = S1[j + 1] + IGST[j + 1];
        }
        
        // Вычисляем T500
        double t500 = RSKVN(PNUM_T500, S1.data(), GI0_DEFAULT);
        
        #pragma omp critical
        {
            *out << std::fixed << std::setprecision(6) << grid[i].first << "\t" << grid[i].second << "\t"
                 << std::fixed << std::setprecision(2) << t500 << std::endl;
            
            if ((i + 1) % 100 == 0 || i + 1 == grid.size()) {
                std::cerr << "Processed " << i + 1 << "/" << grid.size() << " points" << std::endl;
            }
        }
    }
    
    if (fout.is_open()) fout.close();
    
    return 0;
}