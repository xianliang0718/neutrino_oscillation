#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

const float EARTH_RADIUS = 6371.0;

bool firstRun = true;
float degrees_to_decimal(int degrees, int minutes, int seconds)
{
    float deg;
    float min;
    float sec;
    float decimal;

    deg = degrees;
    min = minutes / 60.0;
    sec = seconds / 3600.0;

    decimal = fabs(sec) + fabs(min) + fabs(deg);

    /* 处理负值 */
    if (deg < 0)
    {
        decimal = -decimal;
    }
    else if (deg == 0)
    {
        if (min < 0)
        {
            decimal = -decimal;
        }
        else if (min == 0)
        {
            if (sec < 0)
            {
                decimal = -decimal;
            }
        }
    }

    return (decimal);
}

void calculateDetectorPosition(float latitude, float longitude, float depth, float detectorRadius, float detectorPos[3])
{
    // 将经纬度转换为弧度
    float latRad = latitude * M_PI / 180.0;
    float lonRad = longitude * M_PI / 180.0;

    // 计算探测器中心到地心的距离
    float r = EARTH_RADIUS - depth;

    // 计算直角坐标
    detectorPos[0] = r * cos(latRad) * cos(lonRad);
    detectorPos[1] = r * cos(latRad) * sin(lonRad);
    detectorPos[2] = r * sin(latRad);
}

bool willHitDetector(ROOT::RVec<float> x0, ROOT::RVec<float> p0, float latitude, float longitude, float depth, float R)
{
    // 1. 计算探测器中心位置
    float detectorPos[3];
    calculateDetectorPosition(latitude, longitude, depth, R, detectorPos);

    // 2. 计算粒子运动方向的单位向量
    float pMag = sqrt(p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]);
    if (pMag < 1e-10)
        return false; // 动量过小，视为静止

    float v[3] = {p0[0] / pMag, p0[1] / pMag, p0[2] / pMag};

    for (int i = 0; i < 3; i++)
        x0[i] /= 1e6;


    // 3. 计算向量：粒子位置到探测器中心的向量
    float x0_to_detector[3] = {
        detectorPos[0] - x0[0],
        detectorPos[1] - x0[1],
        detectorPos[2] - x0[2]};

    // 4. 计算最近点参数和距离
    float t = v[0] * x0_to_detector[0] + v[1] * x0_to_detector[1] + v[2] * x0_to_detector[2];

    float distance_Sq = x0_to_detector[0] * x0_to_detector[0] + x0_to_detector[1] * x0_to_detector[1] + x0_to_detector[2] * x0_to_detector[2] - t * t;

    // 5. 判断最近距离是否小于探测器半径
    return distance_Sq <= R * R;
}

void Pick2()
{
    // 定义要提取的 ReactionChain
    std::vector<std::pair<int, std::string>> chains = {
        {12, "ReactionChain(12);363"},
        {12, "ReactionChain(-12);361"},
        {14, "ReactionChain(14);367"},
        {14, "ReactionChain(-14);368"},
    };

    // 创建一个 ROOT 文件用于存储输出
    TFile outputFile("neutrino_data2.root", firstRun ? "RECREATE" : "UPDATE");

    // 创建一个新的 TTree
    TTree outputTree("outputTree", "Filtered Data Tree");

    // 定义需要存储的变量
    int cnt = 0;
    int total = 0;
    int pdgID;
    float Ek0;
    float Zenith;
    float Weight;
    float Height;
    float P0;
    int ilat_deg = 21;
    int ilat_min = 35;
    int ilat_sec = 0;
    int ilon_deg = 112;
    int ilon_min = 30;
    int ilon_sec = 0;
    float latitude, longitude;
    float depth = 0.7;
    float R = 200; // 单位统一km
    ROOT::RVec<float> x0(3);

    // 分配空间给写入树的变量
    outputTree.Branch("pdgID", &pdgID);
    outputTree.Branch("Ek0", &Ek0);
    outputTree.Branch("P0", &P0);
    outputTree.Branch("Zenith", &Zenith);
    outputTree.Branch("Weight", &Weight);
    outputTree.Branch("Height", &Height);

    // 打开一个 txt 文件用于输出
    std::ofstream outputTXT("output1.txt");
    if (!outputTXT.is_open())
    {
        std::cerr << "Failed to open output.txt" << std::endl;
        return;
    }

    for (const auto &chain : chains)
    {
        int currentPdgID = chain.first;
        std::string chainName = "G4Run0/" + chain.second;
        latitude = degrees_to_decimal(ilat_deg, ilat_min, ilat_sec);
        longitude = degrees_to_decimal(ilon_deg, ilon_min, ilon_sec);
        // 打开源 ROOT 文件
        ROOT::RDataFrame d(chainName.c_str(), "/home/xianliang/download/musairs/MusAirS_p_AMS_wo_pipm.root");

        d.Foreach([&cnt, &pdgID, &Ek0, &Zenith, &outputTree, &outputTXT, &Weight, &x0, &Height, &P0, &latitude, &longitude, &R, &depth, &total](int pdgIDValue, float Ek0Value, float ZenithValue, float WeightValue, const ROOT::RVec<float> &xI, const ROOT::RVec<float> PI)
                  {
                    total++;
        bool hit = willHitDetector(xI,PI,latitude,longitude,depth,R);
        if ((pdgIDValue == 12 || pdgIDValue == -12 || pdgIDValue == 14 || pdgIDValue == -14) && hit) {
                cnt++;
                pdgID = pdgIDValue;
                Ek0 = Ek0Value/1000.0;//转化成GeV
                Zenith = ZenithValue;
                Weight = WeightValue;
                for (int i = 0; i < 3; ++i)        x0[i] = xI[i]/1000000.0; // 转换为地心坐标，单位为 km
                Height = sqrt(x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2]) - EARTH_RADIUS;
                outputTree.Fill();
                // 将数据写入 txt 文件
                outputTXT  << pdgID << " " << Ek0 << " " << Zenith << " " << Weight << " " << Height << std::endl;
        } }, {"PDGID", "Ek0", "Zenith", "Weight", "x0", "p0"});
    }
    // 关闭文件和树
    std::cout << "Toatl data = " << total << " Extracted data for all chains, total count: " << cnt << std::endl;
    outputTXT.close();
    outputTree.Write();
    outputFile.Close();
}
