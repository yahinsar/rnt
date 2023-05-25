#define _USE_MATH_DEFINES

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <bitset>
#include <sstream>
#include <deque>
#include <cmath>
#include <set>

using namespace std;

string findInStr(string const& str, int n) {
    if (str.length() < n) {
        return str;
    }
    return str.substr(0, n);
}

void helpFunc() {
    cout << "\n\n/f:<имя_файла> - имя файла с входной последовательностью";
    cout << "\n\n/h – информация о допустимых параметрах командной строки программы.\n";
}


double calculateExpValue(vector <double> allElems) {
    double elemsCount = allElems.size();
    double elemsSum = 0;
    for (int i = 0; i < allElems.size(); i++) {
        elemsSum += allElems[i];
    }
    return elemsSum / elemsCount;
}

double calculateStandartDeviation(vector <double> allElems, double expValue) {
    double elemsCount = allElems.size();
    double elemsSum = 0;
    for (int i = 0; i < allElems.size(); i++) {
        elemsSum = elemsSum + (allElems[i] - expValue) * (allElems[i] - expValue);
    }
    return sqrt(elemsSum / elemsCount);
}

double calculateExpValueError(double expValue) {
    return abs(expValue - 0.5) / expValue;
}

double calculateStandartDeviationError(double standartDeviation) {
    double sqrt1_12 = sqrt(1.0 / 12);
    return abs(standartDeviation - sqrt1_12) / standartDeviation;
}

bool checkXiSquareCriterion(vector <double> allElems) {
    double maxElem = -1.0, minElem = 1.0;
    int elemsCount = allElems.size();
    for (int i = 0; i < elemsCount; i++)
    {
        double checkElem = allElems[i];
        if (checkElem > maxElem)
            maxElem = checkElem;
        if (checkElem < minElem)
            minElem = checkElem;
    }
    vector <pair <double, double>> allIntervals;
    int intervalsCount = ceil(1 + 1.4 * log(elemsCount));
    double IntS = (round(((maxElem - minElem) / intervalsCount) * 10000)) / 10000;
    for (int i = 0; i < intervalsCount; i++) {
        double a_i = minElem + IntS * i;
        double b_i = minElem + IntS * (i + 1);
        allIntervals.push_back(make_pair(a_i, b_i));
    }
    double expectedNumOfHits = round((elemsCount / intervalsCount) * 1000) / 1000;
    vector <int> realDistribution(intervalsCount, 0);
    for (int i = 0; i < elemsCount; i++)
        for (int j = 0; j < intervalsCount; j++)
            if (allIntervals[j].first <= allElems[i] && allElems[i] < allIntervals[j].second) {
                realDistribution[j]++;
                break;
            }
    double XiSquareDistribution = 0;
    for (int i = 0; i < intervalsCount; i++)
        XiSquareDistribution = XiSquareDistribution + ((realDistribution[i] - expectedNumOfHits) * (realDistribution[i] - expectedNumOfHits)) / expectedNumOfHits;
    XiSquareDistribution = round(XiSquareDistribution * 10000) / 10000;
    double criticalImportance = 22.36203;
    cout << "\n > Количество интервалов: " << intervalsCount;
    int k = intervalsCount - 1;
    cout << "\n > Ожидаемое число попаданий в интервалы: " << expectedNumOfHits;
    cout << "\n > Число попаданий в каждый из интервалов: \n";
    for (int i = 0; i < intervalsCount; i++)
        cout << "    " << i + 1 << ": " << realDistribution[i] << endl;
    cout << "\n > Критическое значение хи-квадрат для " << k << " степеней свобод: " << criticalImportance;
    cout << "\n > Значение критерия хи-квадрат с " << k << " степенью свободы: " << XiSquareDistribution << "\n";
    if (0 < XiSquareDistribution && XiSquareDistribution < criticalImportance) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }
}

bool checkSeriesCriterion(vector <double> allElems) {
    int d = 4;
    vector <vector <int>> nPairs(d);
    for (int i = 0; i < d; i++)
        nPairs[i].resize(d, 0);
    int k = d * d;
    vector <pair <int, int>> pairCategories;
    int elemsCount = allElems.size();
    if (elemsCount % 2 != 0)
        elemsCount--;
    for (int j = 0; j < elemsCount; j+=2)
    {
        double leftI = floor(allElems[j] * d);
        double rightI = floor(allElems[j + 1] * d);
        pairCategories.push_back(make_pair(leftI, rightI));
    }
    double criticalImportance = 26.29623;
    double expectedNumOfHits = round((elemsCount / (2 * k)) * 1000) / 1000;
    for (int i = 0; i < pairCategories.size(); i++) {
        nPairs[pairCategories[i].first][pairCategories[i].second] += 1;
    }
    double XiSquareDistribution = 0;
    for (int i = 0; i < nPairs.size(); i++)
        for (int x = 0; x < nPairs[i].size(); x++)
            XiSquareDistribution += ((nPairs[i][x] - expectedNumOfHits) * (nPairs[i][x] - expectedNumOfHits)) / expectedNumOfHits;
    XiSquareDistribution = round(XiSquareDistribution * 1000) / 1000;
    cout << "\n > Параметр d: " << d;
    cout << "\n > Ожидаемое количество чисел в каждой категории: " << expectedNumOfHits;
    cout << "\n > Количество чисел в каждой категории: \n";
    int chCount = 0;
    for (int i = 0; i < nPairs.size(); i++) {
        for (int x = 0; x < nPairs[i].size(); x++) {
            chCount++;
            cout << "   " << chCount << ": " << nPairs[i][x] << "\n";
        }
        chCount++;
    }
    cout << "\n > Критическое значение хи-квадрат для " << k << " степеней свобод: " << criticalImportance;
    cout << "\n > Значение критерия хи-квадрат с " << k << " степенью свободы: " << XiSquareDistribution << "\n";
    if (0 < XiSquareDistribution && XiSquareDistribution < criticalImportance) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }
}

double powDouble(double c, int stepen) {
    if (stepen == 0)
        return 1.0;
    if (stepen == 1)
        return c;
    double res = c;
    for (int i = 2; i <= stepen; i++)
        res = res * c;
    return res;
}

vector < double > InterVector = { 3.84146 , 5.99146, 7.81473, 9.48773, 11.07050, 12.59159, 14.06714, 15.50731, 16.91898, 18.30704, 19.67514, 21.02607, 22.36203, 23.68479, 24.99579, 26.29623, 27.58711, 28.86930, 30.14353, 31.41043, 32.67057, 33.92444, 35.17246, 36.41503, 37.65248, 38.88514, 40.11327, 41.33714, 42.55697, 43.77297 };

bool checkIntervalsCriterion(vector <double> allElems) {
    int t = 10, n = 1000;
    srand((unsigned int)time(0));
    double alpha = (double)(rand()) / RAND_MAX, beta = (double)(rand()) / RAND_MAX;
    while (0.2 >= abs(alpha - beta) || abs(alpha - beta) >= 0.4) {
        alpha = (double)(rand()) / RAND_MAX;
        beta = (double)(rand()) / RAND_MAX;
    }
    if (beta < alpha) {
        double dopD = alpha;
        alpha = beta;
        beta = dopD;
    }
    alpha = round(alpha * 1000) / 1000;
    beta = round(beta * 1000) / 1000;
    double p = beta - alpha;
    int s = 0;
    vector <int> count, resCount;
    double XiSquareDistribution = 100000000000.0;
    double criticalImportance = 100000000000.0;
    int resT = 100000000000;
    int intervalsCount = 100000000000;
    vector <double> resPs;
    s = 0;
    count.resize(t + 1, 0);
    int r = 0;
    for (int j = 0; j < allElems.size(); j++) {
        if (alpha <= allElems[j] && allElems[j] < beta) {
            if (r >= t)
                count[t]++;
            else
                count[r]++;
            s++;
            if (s < n)
                r = 0;
            else
                break;
        }
        else
            r++;
    }
    vector <double> probs;
    for (int j = 0; j < t; j++)
        probs.push_back(round((p * (powDouble(1.0 - p, j))) * 1000) / 1000);
    double XiSch = 0.0;
    probs.push_back(powDouble(1.0 - p, t));
    double criticalImpDop = InterVector[t];
    int qN = 0;
    for (int j = 0; j < count.size(); j++)
        qN += count[j];
    for (int j = 0; j < probs.size(); j++) {
        double est = qN * probs[j];
        if (est == 0.0)
            break;
        XiSch += ((count[j] - est) * (count[j] - est)) / est;
    }
    if (XiSch < XiSquareDistribution) {
        XiSquareDistribution = XiSch;
        criticalImportance = criticalImpDop;
        resT = t;
        intervalsCount = n;
        resPs = probs;
        resCount = count;
    }
    int k = resT + 1;
    cout << "\n > Максимальная длина интервала (t): " << resT;
    cout << "\n > Количество интервалов (n): " << intervalsCount;
    cout << "\n > Границы: alpha = " << alpha << " beta = " << beta;
    cout << "\n > Пересчитанные вероятности Pr и Pt: \n";
    cout << "   ";
    for (int j = 0; j < resPs.size(); j++)
        cout << resPs[j] << " ";
    cout << "\n > Подсчитанные значения интервалов длиной 0, 1, ..., t - 1 и >= t (t = " << resT << "): \n";
    cout << "   ";
    for (int j = 0; j < resCount.size(); j++)
        cout << resCount[j] << " ";
    cout << "\n > Критическое значение хи-квадрат для " << k << " степеней свобод: " << criticalImportance;
    cout << "\n > Значение критерия хи-квадрат с " << k << " степенью свободы: " << round(XiSquareDistribution * 1000) / 1000 << "\n";
    if (0 < XiSquareDistribution && XiSquareDistribution < criticalImportance) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }
}

bool checkPartitionCriterion(vector <double> allElems) {
    int d = 8;
    int k = 5;
    vector <int> resFives(6, 0);
    int elemsCount = allElems.size();
    while (elemsCount % 5 != 0)
        elemsCount = elemsCount - 1;
    int sGr = elemsCount / k;
    for (int i = 0; i < sGr; i++) {
        set <int> dopGr;
        for (int j = i * k; j < (i + 1) * k; j++)
            dopGr.insert(floor(allElems[j] * d));
        resFives[dopGr.size()] += 1;
    }
    double criticalImportance = 9.48773;
    vector <double> fivePs = { 0, 0.0002, 0.0256, 0.2563, 0.5127, 0.2051 };
    double XiSquareDistribution = 0;
    for (int i = 1; i < resFives.size(); i++)
        XiSquareDistribution = XiSquareDistribution + ((resFives[i] - fivePs[i] * sGr) * (resFives[i] - fivePs[i] * sGr)) / (fivePs[i] * sGr);
    XiSquareDistribution = round(XiSquareDistribution * 1000) / 1000;
    cout << "\n > Параметр d: " << d;
    cout << "\n > Параметр k: " << k;
    for (int i = 1; i < resFives.size(); i++) {
        cout << "   Количество пятерок, где " << i << " различных значений: ";
        cout << resFives[i] << "\n";
    }
    cout << "\n > Теоретические вероятности для каждой из пятерок: ";
    for (int i = 1; i < fivePs.size(); i++)
        cout << "   " << i << ": " << fivePs[i] << "\n";
    cout << "\n > Критическое значение хи-квадрат для " << k - 1 << " степеней свобод: " << criticalImportance;
    cout << "\n > Значение критерия хи-квадрат с " << k - 1 << " степенью свободы: " << XiSquareDistribution << "\n";
    if (0 < XiSquareDistribution && XiSquareDistribution < criticalImportance) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }
}

int dopPermutationFunction(vector <double> permVec) {
    int r = 4;
    int fRes = 0;
    while (r > 0) {
        int s = 0;
        double maxV = -10.0;
        for (int i = 0; i < permVec.size(); i++)
            if (maxV < permVec[i]) {
                maxV = permVec[i];
                s = i;
            }
        s++;
        fRes = r * fRes + s - 1;
        double dopV = permVec[r - 1];
        permVec[r - 1] = permVec[s - 1];
        permVec[s - 1] = dopV;
        permVec.resize(permVec.size() - 1);
        r--;
    }
    return fRes;
}

bool checkPermutationsCriterion(vector <double> allElems) {
    int t = 4;
    int myT = 1;
    for (int i = 1; i <= t; i++)
        myT *= i;
    int sizeinVals = allElems.size();
    while (sizeinVals % 4 != 0)
        sizeinVals = sizeinVals - 1;
    double criticalImportance = 36.41503;
    int nGr = sizeinVals / t;
    vector <int> resCategories(myT, 0);
    for (int i = 0; i < nGr; i++) {
        vector <double> cGr;
        set <double> dopGrc;
        for (int j = i * t; j < (i + 1) * t; j++) {
            cGr.push_back(allElems[j]);
            dopGrc.insert(allElems[j]);
        }
        if (dopGrc.size() != t)
            continue;
        else
            resCategories[dopPermutationFunction(cGr)] += 1;
    }
    double expectedP = round((1.0 / myT) * 1000) / 1000;
    double XiSquareDistribution = 0.0;
    double expectedCategoriesValue = round((expectedP * nGr) * 1000) / 1000;
    cout << "\n > Параметр t: " << t;
    cout << "\n > Теоретическое количество попаданий в каждую категорию: " << expectedCategoriesValue;
    cout << "\n > Вычисленное распределение по частотам: ";
    for (int i = 0; i < resCategories.size(); i++)
        cout << resCategories[i] << " ";
    cout << "\n > Теоретическое значение вероятности для каждой категории: " << expectedP;
    for (int i = 0; i < myT; i++)
        XiSquareDistribution = XiSquareDistribution + ((resCategories[i] - expectedCategoriesValue) * (resCategories[i] - expectedCategoriesValue)) / expectedCategoriesValue;
    XiSquareDistribution = round(XiSquareDistribution * 1000) / 1000;
    cout << "\n > Критическое значение хи-квадрат для " << myT << " степеней свобод: " << criticalImportance;
    cout << "\n > Значение критерия хи-квадрат с " << myT << " степенью свободы: " << XiSquareDistribution << "\n";
    if (0 < XiSquareDistribution && XiSquareDistribution < criticalImportance) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }
}

bool checkMonotonicityCriterion(vector <double> allElems) {
    vector <pair <int, int>> categor;
    vector <double> ojid;
    double XiSquareDistribution = 0.0;
    int posl = 1, i = 1;
    while (i < allElems.size() - 1) {
        if (allElems[i] < allElems[i + 1])
            posl += 1;
        else {
            int num = -1;
            if (categor.size() != 0) {
                for (int j = 0; j < categor.size(); j++) {
                    if (categor[j].first == posl) {
                        num = j;
                        break;
                    }
                }
                if (num == -1)
                    categor.push_back(make_pair(posl, 1));
                else
                    categor[num].second += 1;
            }
            else
                categor.push_back(make_pair(posl, 1));
            posl = 1;
            i = i + 1;
        }
        i += 1;
    }
    int sumS = 0;
    for (int j = 0; j < categor.size(); j++)
        sumS = sumS + categor[j].second;
    int num_categs = categor.size();
    double criticalImportance = InterVector[num_categs - 1];
    int maximal = -1;
    for (int j = 0; j < categor.size(); j++)
        if (categor[j].first > maximal)
            maximal = categor[j].first;
    long long fact = 1;
    for (int j = 1; j < maximal + 1; j++) {
        fact *= j;
        ojid.push_back(1.0 / fact - 1.0 / (fact * (j + 1)));
    }
    for (int j = 0; j < num_categs; j++)
        ojid[j] = (sumS * ojid[j]);
    for (int j = 0; j < num_categs; j++)
        XiSquareDistribution = XiSquareDistribution + round((((categor[j].second - ojid[categor[j].first - 1]) * (categor[j].second - ojid[categor[j].first - 1])) / ojid[categor[j].first - 1]) * 1000) / 1000;
    cout << "\n > Теоретические значения длин серий: ";
    for (int j = 0; j < ojid.size(); j++)
        cout << round(ojid[j]) << " ";
    cout << "\n > Вычисленные значения длин серий: ";
    set <pair <int, int>> need;
    for (int j = 0; j < categor.size(); j++)
        need.insert(categor[j]);
    for (pair <int, int> a : need)
        cout << a.first << ": " << a.second << " ";
    cout << "\n > Критическое значение хи-квадрат для " << num_categs << " степеней свобод: " << criticalImportance;
    cout << "\n > Значение критерия хи-квадрат с " << num_categs << " степенью свободы: " << XiSquareDistribution << "\n";
    if (0 < XiSquareDistribution && XiSquareDistribution < criticalImportance) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }

}

vector <pair <int, double>> ConflictCriterionPercentPoints(int m, int n) {
    vector <double> auxiliaryTableT = { 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 1.0 };
    vector <double> auxiliaryA(n + 1, 0.0);
    vector <pair <int, double>> conflictsResP;
    auxiliaryA[1] = 1.0;
    int j0 = 1;
    int j1 = 1;
    for (int i = 0; i < n - 1; i++) {
        j1++;
        for (int j = j1; j > j0 - 1; j--) {
            double jm = j / (m * 1.0);
            auxiliaryA[j] = jm * auxiliaryA[j] + (1.0 + 1.0 / m - jm) * auxiliaryA[j - 1];
            if (auxiliaryA[j] < 1e-20) {
                auxiliaryA[j] = 0.0;
                if (j == j1) {
                    j1--;
                    continue;
                }
                if (j == j0)
                    j0++;
            }
        }
    }
    int t = 0;
    int j = j0 - 1;
    double p = 0.0;
    while (t != auxiliaryTableT.size() - 1) {
        while (p <= auxiliaryTableT[t]) {
            j++;
            p += auxiliaryA[j];
        }
        conflictsResP.push_back(make_pair(n - j - 1, round((1 - p) * 1000) / 1000));
        t++;
    }
    reverse(conflictsResP.begin(), conflictsResP.end());
    return conflictsResP;

}

bool checkConflictCriterion(vector <double> allElems) {
    srand((unsigned int)time(0));
    double eps = 1e-20;
    vector <vector <pair<int, double>>> appropriatePercentPoints;
    vector <vector <int>> appropriatePointsTable;
    for (int bigSize = 8; bigSize < 21; bigSize++) {
        int nParameter = allElems.size() / bigSize;
        for (int dParameter = 2; dParameter < 9; dParameter++) {
            int j = 0;
            vector <int> normSeq;
            while (j < allElems.size()) {
                if (allElems[j] == 1.0)
                    allElems[j] = 0.965;
                normSeq.push_back((floor(allElems[j] * dParameter)));
                j++;
            }
            set <vector <int>> words;
            int nDimension = 0;
            for (int jI = 0; jI < nParameter; jI++) {
                vector <int> dopSliced;
                for (int i = jI * bigSize; i < (jI + 1) * bigSize; i++) {
                    dopSliced.push_back(normSeq[i]);
                }
                auto search = words.find(dopSliced);
                if (search == words.end())
                    words.insert(dopSliced);
                else
                    nDimension += 1;
            }
            for (int mParameter = 16; mParameter < 129; mParameter++) {
                int  mV = nParameter * mParameter;
                vector <pair <int, double>> confs_et_probs = ConflictCriterionPercentPoints(mV, nParameter);
                if (nDimension == 0 || confs_et_probs[0].first == -1 || confs_et_probs[0].first == 0)
                    continue;
                if (confs_et_probs[2].first <= nDimension && nDimension <= confs_et_probs[confs_et_probs.size() - 2].first) {
                    appropriatePercentPoints.push_back(confs_et_probs);
                    vector <int> prom = { nDimension, bigSize, nParameter, dParameter, mParameter, mV };
                    appropriatePointsTable.push_back(prom);
                }
            }
        }
    }
    int appropriatePercentPointsCount = appropriatePercentPoints.size();
    set <int> setRand;
    for (int j = 0; j < 5; j++) {

        int randInt = (double)(rand()) / RAND_MAX * appropriatePercentPointsCount - 1;
        auto search1 = setRand.find(randInt);
        if (search1 == setRand.end())
            setRand.insert(randInt);
        else {
            while (search1 == setRand.end()) {
                randInt = (double)(rand()) / RAND_MAX * appropriatePercentPointsCount - 1;
                search1 = setRand.find(randInt);
            }
        }
    }
    cout << "\n > Число параметров, при которых представленная последовательность удовлетворяет критерию конфликтов: " << appropriatePercentPointsCount;
    cout << "\n > Приведем случай, которая параметры подходят:";
    cout << "\n     Размерность вектора Vj: " << appropriatePointsTable[0][1];
    cout << "\n     Количество векторов: " << appropriatePointsTable[0][2];
    cout << "\n     Значение m: " << appropriatePointsTable[0][5];
    cout << "\n     Множитель для m: " << appropriatePointsTable[0][4];
    cout << "\n     Параметр нормирования d: " << appropriatePointsTable[0][3];
    cout << "\n     Количество возникших конфликтов: " << appropriatePointsTable[0][0];
    cout << "\n     Таблица процентных точек: ";
    for (int k = 0; k < appropriatePercentPoints[0].size(); k++) {
        cout << "(" << appropriatePercentPoints[0][k].first << ", " << appropriatePercentPoints[0][k].second << ") ";
    }
    cout << "\n";
    if (appropriatePercentPointsCount > 0) {
        cout << "\n[+++] Критерий пройден\n";
        return true;
    }
    else {
        cout << "\n[---] Критерий не пройден\n";
        return false;
    }
}

vector <bool> checkCriterions(string inputFileName, bool& errorMessage) {
    vector <bool> criteriaResults(7);
    ifstream inFile(inputFileName);
    string lineFile;
    vector <double> allElems;
    if (!inFile.is_open()) {
        cout << "Ошибка открытия файла " << inputFileName << endl;
        errorMessage = true;
        return criteriaResults;
    }
    string getlineElem;
    getline(inFile, getlineElem);
    stringstream ss(getlineElem);
    while (getline(ss, getlineElem, ',')) {
        if (getlineElem.size() > 1)
            getlineElem[1] = ',';
        double intElem = stod(getlineElem);
        allElems.push_back(intElem);
    }
    inFile.close();

    double expValue = calculateExpValue(allElems);
    cout << "Мат. ожидание: " << round(expValue * 1000) / 1000 << "\n";
    double standartDeviation = calculateStandartDeviation(allElems, expValue);
    cout << "Среднекв. отклонение: " << round(standartDeviation * 1000) / 1000 << "\n";
    double expValueError = calculateExpValueError(expValue);
    cout << "Погрешность мат. ожидания: " << round(expValueError * 1000) / 1000 << "\n";
    double standartDeviationError = calculateStandartDeviationError(standartDeviation);
    cout << "Погрешность среднекв. отклонения: " << round(standartDeviationError * 1000) / 1000 << "\n";

    bool res = false;
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий Хи-квадрат \n\n";
    criteriaResults[0] = checkXiSquareCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий серий \n\n";
    criteriaResults[1] = checkSeriesCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий интервалов \n\n";
    criteriaResults[2] = checkIntervalsCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий разбиений \n\n";
    criteriaResults[3] = checkPartitionCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий перестановок \n\n";
    criteriaResults[4] = checkPermutationsCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий монотонности \n\n";
    criteriaResults[5] = checkMonotonicityCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";
    cout << "\n Критерий конфликтов \n\n";
    criteriaResults[6] = checkConflictCriterion(allElems);
    cout << "\n ----------------------------------------------------------------------- \n";

    return criteriaResults;
}

void coutCriterionsResult(vector <bool> criteriaResults) {
    vector <string> criteriaNames = { "Критерий хи-квадрат", "Критерий серий",  "Критерий интервалов", "Критерий разбиений",  "Критерий перестановок",  "Критерий монотонности",  "Критерий конфликтов"};
    cout << "\nРезультат тестирования статистических свойств последовательности ПСЧ: \n\n";
    for (int i = 0; i < criteriaResults.size(); i++)
    {
        cout << criteriaNames[i] << ": ";
        if (criteriaResults[i])
            cout << "+++" << endl;
        else
            cout << "---" << endl;
    }
    return;
}

int main(int argc, char* argv[])
{
    /*
            Программа для тестирования статистических свойств последовательности ПСЧ
    */
    setlocale(LC_ALL, "Rus");
    string inputFileName = "rntInput.dat";
    for (int i = 0; argv[i]; i++)
    {
        string checkStr = string(argv[i]);
        if (findInStr(checkStr, 2) == "/h") {
            helpFunc();
            return 0;
        }
        if (checkStr.length() > 2) {
            string ifStr = findInStr(checkStr, 3);
            string subStr = checkStr.substr(3, checkStr.length());
            if (ifStr == "/f:") {
                inputFileName = subStr;
            }
        }
    }
    bool errorMessage = false;
    vector <bool> resString = checkCriterions(inputFileName, errorMessage);
    if (errorMessage)
        return 0;

    cout << "Input File: " << inputFileName << endl;

    coutCriterionsResult(resString);
    return 0;
}
