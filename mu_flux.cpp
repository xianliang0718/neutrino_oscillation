#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <cmath>
#include "TH2D.h"
#include "TFile.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <nuSQUIDS/nuSQuIDS.h>

using namespace nusquids;
using namespace std;

// 定义一个结构体来存储数据点的x和y值
struct DataPoint
{
    double x;
    double y;
};

// 从CSV文件中读取数据点并存储在vector中
vector<DataPoint> readCSV(const string &fileName)
{
    vector<DataPoint> data;
    ifstream file(fileName);
    if (!file.is_open())
    {
        cerr << "Error: file not found " << fileName << endl;
        return data;
    }
    string line;
    while (getline(file, line))
    {
        stringstream ss(line);
        string x_str, y_str;
        getline(ss, x_str, ',');
        getline(ss, y_str, ',');
        DataPoint point;
        point.x = stod(x_str);
        point.y = stod(y_str);
        data.push_back(point);
    }
    return data;
}

// 根据读取的数据和自定义的bin中心计算相对值，并存储在map中
map<double, double> getZenithRelVal(const vector<DataPoint> &data, const vector<double> &binCenters)
{
    map<double, double> relValCosZenith;
    if (data.empty())
    {
        cerr << "Error: no data read from file" << endl;
        return relValCosZenith;
    }
    const int min_points = gsl_interp_type_min_size(gsl_interp_cspline);
    if (data.size() < min_points)
    {
        cerr << "Error: Need at least " << min_points
             << " data points. Only " << data.size() << " found." << endl;
        return relValCosZenith;
    }
    for (size_t i = 1; i < data.size(); ++i)
    {
        if (data[i].x <= data[i - 1].x)
        {
            cerr << "Error: x values must be strictly increasing." << endl;
            return relValCosZenith;
        }
    }

    // 初始化GSL插值加速器和样条插值器
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, data.size());

    // 将数据点的x和y值分别存储在两个vector中
    vector<double> xVals(data.size()), yVals(data.size());
    for (size_t i = 0; i < data.size(); i++)
    {
        xVals[i] = data[i].x;
        yVals[i] = data[i].y;
    }

    // 初始化样条插值器
    gsl_spline_init(spline, xVals.data(), yVals.data(), data.size());

    // 遍历自定义的bin中心
    for (double cz : binCenters)
    {
        if (cz < xVals.front() || cz > xVals.back())
        {
            cerr << "Warning: bin center " << cz << " is out of data range." << endl;
            continue;
        }
        double flux = gsl_spline_eval(spline, cz, acc);
        relValCosZenith[cz] = flux;
    }

    // 找到relValCosZenith中的最大值
    double maxVal = 0.0;
    for (const auto &pair : relValCosZenith)
    {
        if (pair.second > maxVal)
        {
            maxVal = pair.second;
        }
    }

    // 用最大值归一化relValCosZenith
    for (auto &pair : relValCosZenith)
    {
        pair.second /= maxVal;
    }

    // 释放GSL插值器和加速器
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return relValCosZenith;
}

// 根据读取的数据和自定义的bin中心计算能量光谱，并存储在map中
map<double, double> getEnergySpectrum(const vector<DataPoint> &data, const vector<double> &binCenters) {
    map<double, double> relValEnergy; // 用于存储归一化的能量光谱

    // 检查数据是否为空
    if (data.empty()) {
        cerr << "Error: no data read from file" << endl;
        return relValEnergy;
    }

    // 获取进行样条插值所需的最小数据点数量
    const int min_points = gsl_interp_type_min_size(gsl_interp_cspline);
    if (data.size() < min_points) {
        cerr << "Error: Need at least " << min_points
             << " data points. Only " << data.size() << " found." << endl;
        return relValEnergy;
    }

    // 检查数据点的x值是否严格递增
    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i].x <= data[i - 1].x) {
            cerr << "Error: x values must be strictly increasing." << endl;
            return relValEnergy;
        }
    }

    // 初始化GSL插值加速器和样条插值器
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, data.size());

    // 将数据点的x和y值分别存储在两个vector中
    vector<double> xVals(data.size()), yVals(data.size());
    for (size_t i = 0; i < data.size(); i++) {
        xVals[i] = data[i].x;
        yVals[i] = data[i].y;
    }

    // 初始化样条插值器
    gsl_spline_init(spline, xVals.data(), yVals.data(), data.size());

    // 遍历自定义的bin中心
    for (double binCenter : binCenters) {
        auto center = log10(binCenter); // 将bin中心转换为对数形式
        if (center < xVals.front() || center > xVals.back()) {
            cerr << "Warning: bin center " << binCenter << " is out of data range." << endl;
            continue; // 如果bin中心超出数据范围，跳过该bin中心
        }
        double flux = gsl_spline_eval(spline, center, acc); // 使用样条插值计算flux
        relValEnergy[binCenter] = flux; // 将计算得到的flux存储在map中
    }

    // 释放GSL插值器和加速器
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return relValEnergy;
}


// 生成中微子通量并进行演化
void generateFlux(double theta32, int theta32Num, double deltaM13, int deltaM13Num)
{
    // 定义单位和常数类
    squids::Const units;
    unsigned int numneu;
    numneu = 3;               // 设置中微子类型数量为3（通常为电子中微子、μ中微子和τ中微子）
    bool interactions = true; // 启用相互作用

    // 设置能量最小值和最大值，以及cos(zenith)的最小值和最大值
    double Emin = 1.e-1 * units.GeV; // 设置能量最小值为0.1 GeV
    double Emax = 1.e2 * units.GeV;  // 设置能量最大值为100 GeV
    double czmin = -1;               // 设置cos(zenith)最小值为-1
    double czmax = 1;                // 设置cos(zenith)最大值为1

    // 构造nuSQUIDSAtm对象
    std::cout << "Begin: constructing nuSQUIDS-Atm object" << std::endl;
    nuSQUIDSAtm<> nus_atm(
        linspace(czmin, czmax, 40), // 生成40个等间距的cos(zenith)值
        logspace(Emin, Emax, 100),  // 生成100个对数间距的能量值
        numneu,                     // 中微子类型数量
        both,                       // 同时考虑中微子和反中微子
        interactions                // 启用相互作用
    );
    std::cout << "End: constructing nuSQUIDS-Atm object" << std::endl;

    // 从CSV文件中读取中微子和反中微子的能量和cos(zenith)光谱
    vector<DataPoint> numu_energy_spectrum = readCSV("/home/smc/software/nuSQuIDS/data/hkkm/real_atn_flux_numu.csv");
    vector<DataPoint> nue_energy_spectrum = readCSV("/home/smc/software/nuSQuIDS/data/hkkm/real_atn_flux_nue.csv");
    vector<DataPoint> numu_cosZenith_spectrum = readCSV("/home/smc/software/nuSQuIDS/data/hkkm/atn_zenith_numu.csv");
    vector<DataPoint> nue_cosZenith_spectrum = readCSV("/home/smc/software/nuSQuIDS/data/hkkm/atn_zenith_nue.csv");

    // 设置混合角、平方质量差和CP相位
    std::cout << "Begin: setting mixing angles." << std::endl;
    nus_atm.Set_MixingAngle(0, 1, 0.563942);
    nus_atm.Set_MixingAngle(0, 2, 0.154085);
    nus_atm.Set_MixingAngle(1, 2, theta32); // 使用输入参数设置θ23混合角
    nus_atm.Set_SquareMassDifference(1, 7.65e-05);
    nus_atm.Set_SquareMassDifference(2, deltaM13); // 使用输入参数设置Δm31^2
    nus_atm.Set_CPPhase(0, 2, 0);
    if (numneu > 3)
    {
        nus_atm.Set_SquareMassDifference(3, -1.);
        nus_atm.Set_MixingAngle(1, 3, 0.160875);
    }
    std::cout << "End: setting mixing angles." << std::endl;

    // 设置积分精度
    nus_atm.Set_rel_error(1.0e-6);             // 设置相对误差
    nus_atm.Set_abs_error(1.0e-6);             // 设置绝对误差
    nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4); // 设置GSL步进方法

    // 获取能量和cos(zenith)范围
    auto e_range = nus_atm.GetERange();
    auto cz_range = nus_atm.GetCosthRange();

    // 将能量和cos(zenith)范围转换为vector
    vector<double> eVec;
    vector<double> czVec;
    for (const auto &e : e_range)
    {
        eVec.push_back(e);
    }
    for (const auto &cz : cz_range)
    {
        czVec.push_back(cz);
    }

    // 计算归一化的能量和cos(zenith)光谱
    map<double, double> numu_energy = getEnergySpectrum(numu_energy_spectrum, eVec);
    map<double, double> nue_energy = getEnergySpectrum(nue_energy_spectrum, eVec);
    map<double, double> numu_cosZenith = getZenithRelVal(numu_cosZenith_spectrum, czVec);
    map<double, double> nue_cosZenith = getZenithRelVal(nue_cosZenith_spectrum, czVec);

    // 设置初始状态
    std::cout << "Begin: setting initial state." << std::endl;
    // 构造初始状态，假设在味道基下设置初始状态
    marray<double, 4> inistate{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, numneu};
    std::fill(inistate.begin(), inistate.end(), 0); // 初始化所有分量为0

    // 遍历所有cos(zenith)、能量、密度状态和中微子类型
    for (int ci = 0; ci < nus_atm.GetNumCos(); ci++)
    {
        for (int ei = 0; ei < nus_atm.GetNumE(); ei++)
        {
            for (int rho = 0; rho < 2; rho++)
            { // rho=0表示中微子，rho=1表示反中微子
                for (int flv = 0; flv < numneu; flv++)
                { // 遍历所有中微子类型
                    if (flv == 0)
                    { // 电子中微子
                        inistate[ci][ei][rho][flv] = nue_energy[e_range[ei]] * nue_cosZenith[cz_range[ci]];
                    }
                    else if (flv == 1)
                    { // μ中微子
                        inistate[ci][ei][rho][flv] = numu_energy[e_range[ei]] * numu_cosZenith[cz_range[ci]];
                    }
                    else
                    {                                     // 其他中微子类型（如τ中微子）
                        inistate[ci][ei][rho][flv] = 0.0; // 暂时设置为0
                    }
                }
            }
        }
    }

    // 在atmSQuIDS对象中设置初始状态
    nus_atm.Set_initial_state(inistate, flavor);
    std::cout << "End: setting initial state." << std::endl;

    // 设置进度条和真空振荡
    nus_atm.Set_ProgressBar(true);         // 启用进度条
    nus_atm.Set_IncludeOscillations(true); // 启用振荡

    // 进行状态演化
    std::cout << "Begin: Evolution" << std::endl;
    nus_atm.EvolveState();
    std::cout << "End: Evolution" << std::endl;

    // 将最终状态保存到文本文件中
    string outputFileName = "atn_" + to_string(theta32Num) + "_" + to_string(deltaM13Num) + ".txt";
    std::ofstream file(outputFileName);
    file << "# log10(E) cos(zenith) E flux_i . . . ." << std::endl;

    // 设置输出的分辨率和范围
    int Nen = 700; // 能量分辨率
    int Ncz = 100; // cos(zenith)分辨率
    double lEmin = log10(Emin);
    double lEmax = log10(Emax);

    // 写入文件
    for (double cz = czmin; cz < czmax; cz += (czmax - czmin) / (double)Ncz)
    {
        for (double lE = lEmin; lE < lEmax; lE += (lEmax - lEmin) / (double)Nen)
        {
            double E = pow(10.0, lE); // 将对数能量转换为实际能量
            file << lE - log10(units.GeV) << " " << cz << " " << E / units.GeV;
            for (int fl = 0; fl < numneu; fl++)
            {
                file << " " << nus_atm.EvalFlavor(fl, cz, E, 0); // 评估中微子在能量和cos(zenith)下的通量
            }
            for (int fl = 0; fl < numneu; fl++)
            {
                file << " " << nus_atm.EvalFlavor(fl, cz, E, 1); // 评估反中微子在能量和cos(zenith)下的通量
            }
            file << std::endl;
        }
        file << std::endl;
    }
}

// 主函数，接收命令行参数
int main(int argc, char *argv[])
{
    // 检查输入参数数量是否足够
    if (argc < 5)
    {
        cerr << "Error: not enough input parameters." << endl;
        cerr << "Usage: " << argv[0] << " <mixing angle> <mixing angle number> <mass difference> <mass difference number>" << endl;
        return 1;
    }

    // 从命令行参数中读取混合角和质量差
    double mixingAngle = stod(argv[1]);    // θ23混合角
    int mixingAngleNum = stoi(argv[2]);    // 混合角编号
    double massDifference = stod(argv[3]); // Δm31^2
    int massDifferenceNum = stoi(argv[4]); // 质量差编号

    // 调用generateFlux函数生成中微子通量并进行演化
    generateFlux(mixingAngle, mixingAngleNum, massDifference, massDifferenceNum);

    return 0;
}
