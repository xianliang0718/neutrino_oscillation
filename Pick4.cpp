#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>

// 手动实现 gsl_vector_bsearch
size_t gsl_vector_bsearch(const gsl_vector *v, double x) {
    size_t low = 0;
    size_t high = v->size;

    while (low < high) {
        size_t mid = low + (high - low) / 2;
        if (gsl_vector_get(v, mid) < x) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

void Pick4() {
    std::ifstream inputTXT("output2.txt");
    if (!inputTXT.is_open()) {
        std::cerr << "无法打开输入文件!" << std::endl;
        return;
    }

    int pdgIDValue;
    double Ek0Value;
    double ZenithValue;
    double WeightValue;
    double Emin = 1.e-3; 
    double Emax = 1.e6; 
    const double GeV_to_eV = 1e9; // 1 GeV = 1e9 eV
    std::map<int, double> e_W_Map; // 用于存储每种 pdgID 对应的总 Weight

    size_t E_num = 200; 
    gsl_vector *enu_marray = gsl_vector_alloc(E_num);

    // 创建对数空间数组
    for (size_t i = 0; i < E_num; ++i) {
        double logE = log(Emin * GeV_to_eV) + (i * (log(Emax * GeV_to_eV) - log(Emin * GeV_to_eV)) / (E_num - 1));
        gsl_vector_set(enu_marray, i, logE);
    }

    // 确保数组已排序
    double *data = gsl_vector_ptr(enu_marray, 0);
    gsl_sort(data, 1, E_num);

    while (inputTXT >> pdgIDValue >> Ek0Value >> ZenithValue >> WeightValue) {
        // 累加每种 pdgID 的 Weight
        double logE = log(Ek0Value * GeV_to_eV);

        // 使用gsl_vector_bsearch查找插入位置
        size_t enu_M_1 = gsl_vector_bsearch(enu_marray, logE);
        size_t enu_M_0 = (enu_M_1 > 0) ? enu_M_1 - 1 : 0; // 防止出界

        // 获取这些位置对应的值
        double enu_m1_value = gsl_vector_get(enu_marray, enu_M_1);
        double enu_m0_value = gsl_vector_get(enu_marray, enu_M_0);

        e_W_Map[pdgIDValue] += WeightValue * (enu_m1_value - enu_m0_value);
    }

    inputTXT.close();

    double totalWeight = 0.0;
    // 输出累加的结果到控制台
    for (const auto &entry : e_W_Map) {
        pdgIDValue = entry.first;
        WeightValue = entry.second;
        std::cout << "pdgID: " << pdgIDValue << ", Total Weight: " << WeightValue << std::endl;
        totalWeight += WeightValue;
    }
    std::cout << "Total Weight: " << totalWeight << std::endl;

    gsl_vector_free(enu_marray);
}

int main() {
    Pick4();
    return 0;
}
