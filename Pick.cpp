#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

bool firstRun = true;

void Pick()
{
    // 定义要提取的 ReactionChain
    std::vector<std::pair<int, std::string>> chains = {
        {12, "ReactionChain(12);34"},
        {12, "ReactionChain(12);33"},
        {-12, "ReactionChain(-12);33"},
        {-12, "ReactionChain(-12);32"},
        {14, "ReactionChain(14);38"},
        {14, "ReactionChain(14);37"},
        {-14, "ReactionChain(-14);38"},
        {-14, "ReactionChain(-14);37"}
    };

    int cnt = 0;

    // 创建一个 ROOT 文件用于存储输出
    TFile outputFile("neutrino_data1.root", firstRun ? "RECREATE" : "UPDATE");

    // 创建一个新的 TTree
    TTree outputTree("outputTree", "Filtered Data Tree");

    // 定义需要存储的变量
    int pdgID;
    ROOT::RVec<float> x0(3);
    ROOT::RVec<float> x(3);
    float Ek0,t,Zenith;

    // 为 TTree 创建分支
    outputTree.Branch("PDGID", &pdgID);
    outputTree.Branch("x0", &x0);
    outputTree.Branch("x", &x);
    outputTree.Branch("Ek0", &Ek0);
    outputTree.Branch("t", &t);
    outputTree.Branch("Zenith", &Zenith);

    // 遍历各个 ReactionChain
    for (const auto &chain : chains)
    {
        int currentPdgID = chain.first;
        std::string chainName = "G4Run0/" + chain.second;

        // 打开源 ROOT 文件
        ROOT::RDataFrame d(chainName.c_str(), "/home/xianliang/download/musairs/MusAirS_p_AMS_wo_pipm.root");

        // 使用 RDataFrame 过滤数据并填充 TTree
        d.Foreach([&cnt, &pdgID, &x0, &x, &Ek0, &t, &outputTree, &Zenith](int pdgIDValue, const ROOT::RVec<float> &xI, const ROOT::RVec<float> &xF, float Ek0Value, float tValue, float ZenithValue)
                  {
            if (pdgIDValue == 12 || pdgIDValue == -12 || pdgIDValue == 14 || pdgIDValue == -14)
            {
                cnt++;
                pdgID = pdgIDValue;
                    for (int i = 0; i < 3; ++i) {
                        x0[i] = xI[i];
                        x[i] = xF[i];
                    }
                Ek0 = Ek0Value;
                t = tValue;
                Zenith = ZenithValue;
                outputTree.Fill();
            } }, {"PDGID", "x0", "x", "Ek0","t", "Zenith"});

        std::cout << "Extracted data for " << chainName << "  cnt=" << cnt << std::endl;
    }

    // 写入 TTree 到文件
    outputTree.Write();
    outputFile.Close();

    // 统计结果
    std::cout << "Total count of matched PDGID: " << cnt << std::endl;
    std::cout << "Data extraction complete!" << std::endl;

    // 直接绘制 PDGID 为 ±12 和 ±14 的数据到直方图中
    // TH1F *histEk0 = new TH1F("histEk0", "Initial Energy Histogram for neutrinos", 10000, 100, 1000000);

    // // 创建 RDataFrame 读取 outputTree 数据
    // ROOT::RDataFrame df("outputTree", "neutrino_data1.root"); // 修改文件名以匹配输出文件

    // // 用 Filter 筛选 PDGID，并填充直方图
    // df.Foreach([&histEk0](int pdgIDValue, float x) {
    //   // 这里是对过滤后的数据进行填充
    //   if(pdgIDValue == 12 || pdgIDValue== -12 || pdgIDValue == 14 || pdgIDValue == -14|| pdgIDValue == 16 || pdgIDValue == -16)
    //          histEk0->Fill(x);
    //          }, {"PDGID", "x"});

    // // 绘制直方图
    // TCanvas *c = new TCanvas("c", "Energy Distribution", 800, 600);
    // c->SetLogx(); // 设置横坐标为 log10 形式
    // histEk0->Draw();
    // c->SaveAs("outputTree_energy_histogram_pdgid.png"); // 保存直方图为图片

    // // 清理
    // delete histEk0;
}
