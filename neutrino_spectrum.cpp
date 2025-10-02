#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>

bool firstRun = true;

void neutrino_spectrum()
{
    std::ifstream inputTXT("/home/xianliang/download/musairs/output1.txt");
    if (!inputTXT)
    {
        std::cerr << "无法打开输入文件" << std::endl;
        return;
    }

    TFile outputFile("neutrino_initial.root", firstRun ? "RECREATE" : "UPDATE");
    TTree outputTree("outputTree", "Initial Data Tree");

    int pdgID;
    double Ek0;
    double Zenith;
    double Weight;

    outputTree.Branch("pdgID", &pdgID);
    outputTree.Branch("Ek0", &Ek0);
    outputTree.Branch("Zenith", &Zenith);
    outputTree.Branch("Weight", &Weight);

    int pdgIDValue;
    double Ek0Value;
    double ZenithValue;
    double WeightValue;
    double HeightValue;
    std::map<int, std::vector<std::pair<double, double>>> e_W_Map;

    int nBins = 20;
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::lowest();
    double Total_Weight = 0;

    while (inputTXT >> pdgIDValue >> Ek0Value >> ZenithValue >> WeightValue >> HeightValue)
    {
        if (Ek0Value > 1e-1 && Ek0Value < 1e3)
        {
            if (Ek0Value < xMin)
                xMin = Ek0Value;
            if (Ek0Value > xMax)
                xMax = Ek0Value;
            e_W_Map[std::abs(pdgIDValue)].emplace_back(Ek0Value, WeightValue*pow(Ek0Value, 2));
            outputTree.Fill();
        }
    }

    std::cout << "xMax=" << xMax << " xMin=" << xMin << std::endl;
    xMax *= 5;
    xMin /= 5;
    inputTXT.close();

    // 对每个 pdgID 的数据按横坐标（能量）排序
    for (auto &entry : e_W_Map)
    {
        std::sort(entry.second.begin(), entry.second.end(), [](const std::pair<double, double> &a, const std::pair<double, double> &b)
                  {
                      return a.first < b.first; // 按能量排序
                  });
    }

    outputTree.Write();
    outputFile.Close();

    std::map<int, TH1D *> Bin_Hist_Map;
    double bins[nBins + 1];
    for (int i = 0; i <= nBins; i++)
    {
        bins[i] = TMath::Power(10, TMath::Log10(xMin) + (TMath::Log10(xMax) - TMath::Log10(xMin)) * i / nBins);
    }

    // 初始化直方图
    for (const auto &entry : e_W_Map)
    {
        int targetPdgID = entry.first;
        // 创建 TH1D 对象，参数分别是名称、标题、bin 的数量、xMin、xMax
        std::string histName = "hist_" + std::to_string(targetPdgID);
        TH1D *hist = new TH1D(histName.c_str(), "Energy Spectrum;Energy (GeV);flux * E_{#nu}^{2}(m^{-1}s^{-1}sr^{-1}GeV)", nBins, bins);
        Bin_Hist_Map[targetPdgID] = hist;
    }

    // 填充数据到直方图
    double Total_Flux = 0;
    for (const auto &entry : e_W_Map)
    {
        const auto &data = entry.second;
        int targetPdgID = std::abs(entry.first);
        int cnt = 0;
        double weight = 0;
        for (const auto &point : data)
        {
            Bin_Hist_Map[targetPdgID]->Fill(point.first, point.second); // TH1D 的 Fill 函数直接填充能量和权重
            cnt++;
            weight += point.second;
        }
        Total_Flux += weight;
        cout << "pdgID: " << targetPdgID << " count: " << cnt << " weight: " << weight << endl;
    }

    std::cout << "Total_Flux: " << Total_Flux << std::endl;
    
    // 创建画布
    TCanvas *canvas = new TCanvas("canvas", "Energy Spectrum", 800, 600);
    canvas->SetGrid();
    canvas->SetLogx(); // 设置x轴为对数刻度

    TLegend *legend = new TLegend(0.1, 0.7, 0.25, 0.9);
    legend->SetHeader("neutrino flavors");

    int colorIndex = 1;                                  // 用于为每个 pdgID 分配不同的颜色
    double yMin = std::numeric_limits<double>::max();    // 初始化 yMin 为最大值
    double yMax = std::numeric_limits<double>::lowest(); // 初始化 yMax 为最小值

    // 填充 yMin 和 yMax
    for (const auto &entry : Bin_Hist_Map)
    {
        const auto &hist = entry.second;
        for (int bin = 1; bin <= hist->GetNbinsX(); bin++)
        { // 使用 GetNbinsX() 获取直方图的总箱数
            double content = hist->GetBinContent(bin);
            if (content > 0)
            { // 确保内容为正
                if (content < yMin)
                    yMin = content;
                if (content > yMax)
                    yMax = content;
            }
        }
    }

    // 调整 yMin 和 yMax 以适应比例尺
    yMin /= 1.1;
    yMax *= 1.1;

    // 绘制直方图
    bool firstPlot = true;
    for (const auto &entry : Bin_Hist_Map)
    {
        int targetPdgID = entry.first;
        TH1D *hist = entry.second;

        hist->SetLineColor(colorIndex);
        hist->SetLineWidth(2);
        hist->SetFillColor(0); // 不填充颜色

        if (firstPlot)
        {
            hist->Draw("HIST"); // 第一个图绘制轴和线
            firstPlot = false;
        }
        else
        {
            hist->Draw("HIST SAME"); // 其余的只绘制线
        }

        // 设置y轴范围
        hist->GetYaxis()->SetRangeUser(yMin, yMax);

        // 添加图例条目
        if (targetPdgID == 12)
            legend->AddEntry(hist, Form("#nu_{e}"), "l");
        else if (targetPdgID == -12)
            legend->AddEntry(hist, Form("#bar{#nu}_{e}"), "l");
        else if (targetPdgID == 14)
            legend->AddEntry(hist, Form("#nu_{#mu}"), "l");
        else if (targetPdgID == -14)
            legend->AddEntry(hist, Form("#bar{#nu}_{#mu}"), "l");
        else if (targetPdgID == 16)
            legend->AddEntry(hist, Form("#nu_{#tau}"), "l");
        else if (targetPdgID == -16)
            legend->AddEntry(hist, Form("#bar{#nu}_{#tau}"), "l");

        // 调整图例的文本大小
        legend->SetTextSize(0.04); // 增加文本大小
        // 设置文本对齐方式（可选）
        legend->SetTextAlign(12); // 左对齐，垂直居中对齐
        // 设置图例框的边框大小和填充颜色（可选）
        legend->SetBorderSize(1);
        legend->SetFillColor(0);

        colorIndex++;
        if (colorIndex > 8)
            colorIndex = 1;
    }

    // 绘制图例
    legend->Draw();

    // 保存画布为图片文件
    canvas->SaveAs("neutrino_initial.png");
}

int main()
{
    neutrino_spectrum();
    return 0;
}
// th1d