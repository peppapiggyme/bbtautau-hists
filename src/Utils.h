#ifndef UTILS_H
#define UTILS_H

#include "Processes.h"
#include "Regions.h"
#include "Variables.h"
#include "Systematics.h"

#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <climits>

using std::string;
using std::vector;
using std::pair;

namespace Utils {

    enum class NameConvention { CxAODReader, WSMaker };

    void histAssign(TH1* h, ProcessInfo* p, RegionInfo* r, VariableInfo* v);
    void histAssignSyst(TH1* h, ProcessInfo *p, const std::string& systname);

    /**
     * @note so far only support CxAODReader style naming
     */
    string histString(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v, const NameConvention nc = NameConvention::CxAODReader);
    
    /**
     * @note so far only support CxAODReader style naming
     */
    string systString(const SystematicInfo* s, const NameConvention nc = NameConvention::CxAODReader);
    
    /**
     * @note so far only support CxAODReader style naming
     */
    string histStringSyst(const ProcessInfo* p, const RegionInfo* r, const VariableInfo* v, const SystematicInfo* s, const NameConvention nc = NameConvention::CxAODReader);
    
    /**
     * @brief shorten the names of nuisance parameters 
     */
    string systStringShort(const string& sSyst, const NameConvention nc = NameConvention::CxAODReader);
    
    /**
     * @brief histogram process name are hhttbb, Hhhbbtautau, ...
     * In WSMaker, the signal type name are SM, 2HDM, ...
     */
    string signalTypeName(const string& sSigName);

    /**
     * @brief the standard directory of background-only fit
     * i.e. the POI (signal strength) was set to constant of 0 in the fit
     */
    inline const string bkgOnlyFitResultString() { return "PlotsAfterGlobalFit/conditionnal_MuIsEqualTo_0/fitResult"; }

    void properties_copy(TH1* h1, TH1* h2);
}

class BinningUtils
{
public:

    template<typename T>
    static std::vector<T> readBinningFromFile(const std::string& fn)
    {
        std::vector<T> binEdges;

        std::ifstream f(fn);
        std::string line;
        
        getline(f, line);
        
        std::istringstream iss(line);
        T x;
        while (iss >> x)
        {
            binEdges.push_back(x); // off by one
        }

        return binEdges;
    }

    // static std::vector<double> intToDoubleBinEdgesForMVA(const std::vector<int>& ii, int nbins=1000)
    // {
    //     std::vector<double> ret(ii.size()+1);
    //     ret.at(0) = 0.;

    //     for (size_t i = 1; i <= ii.size(); ++i)
    //     {
    //         ret.at(i) = 0. + (double)ii[i]/(double)nbins;
    //     }
        
    //     return ret;
    // }

    static std::vector<double> intToDoubleBinEdgesForMVAInverse(const std::vector<int>& ii, int nbins=1000, bool isBDT=false)
    {
        std::vector<double> ret(ii.size());
        if (isBDT) 
        {
            ret.at(0) = -1.;
        } 
        else
        {
            ret.at(0) = 0.;
        }
        ret.at(ret.size()-1) = 1.;

        for (size_t i = 1; i < ret.size()-1; ++i)
        {
            int ibin = ii[ret.size()-1-i] - 1;
            if (isBDT)
            {
                ret.at(i) = -1. + 2. * (double)ibin/(double)nbins;
            }
            else
            {
                ret.at(i) = 0. + (double)ibin/(double)nbins;
            }
        }
        
        return ret;
    }

    static std::vector<double> intToDoubleBinEdgesForMVAInverse_LepHad1090Version(const std::vector<int>& ii)
    {
        std::vector<double> ret(ii.size());

        ret.at(0) = 0.;
        
        ret.at(ret.size()-1) = 1.;

        for (size_t i = 1; i < ret.size()-1; ++i)
        {
            int ibin = ii[ret.size()-1-i] - 1;
            if (ibin <= 990)
                ret.at(i) = 0. + (double)ibin / 1e3;
            else 
                ret.at(i) = 0.99 + (ibin - 990) / 1e4;
        }
        
        return ret;
    }

};

class Tools
{

public:
    static void print(const char* fmt)
    {
        std::cout << fmt;
    }

    static void println(const char* fmt)
    {
        Tools::print(fmt);
        Tools::print("\n");
    }

    template<typename T, typename ... Targs>
    static void print(const char* fmt, T value, Targs ... args)
    {
        for(; *fmt != '\0'; fmt++)
        {
            if (*fmt == '%') {
                std::cout << value;
                Tools::print(fmt+1, args...);
                return;
            }
            std::cout << *fmt;
        }
    }

    template<typename T, typename ... Targs>
    static void println(const char* fmt, T value, Targs ... args)
    {
        Tools::print(fmt, value, args...);
        Tools::print("\n");
    }

    static int getInteger(const std::string& prompt = "Type in an integer: ", 
                          const std::string& reprompt = "It is not integer. Retry. \n")
    {
        while (1)
        {
            std::cout << prompt;
            std::string line;
            if (!getline(std::cin, line)) throw std::domain_error("Failed to get line from cin.");

            std::istringstream iss(line);
            int i; char a;
            if (iss >> i && !(iss >> a))
            {
                return i;
            }
            std::cout << reprompt;
        }
        return INT_MIN;
    }

    static std::string getString(const std::string& prompt = "Type in an string: ", 
                          const std::string& reprompt = "It is not string. Retry. \n")
    {
        while (1)
        {
            std::cout << prompt;
            std::string line;
            if (!getline(std::cin, line)) throw std::domain_error("Failed to get line from cin.");

            std::istringstream iss(line);
            std::string str; char a;
            if (iss >> str && !(iss >> a))
            {
                return str;
            }
            std::cout << reprompt;
        }
        return "";
    }

    /**
     * @brief pretty print vectors
     * @note pass by value so the print uses a copy
     */
    template<typename T>
    static void printVector(const T& v) 
    {
        size_t cnt = v.size();
        for (const auto& x : v)
        {
            std::cout << x;
            if (cnt != 1)
            {
                std::cout << ", ";
            }
            cnt--;
        }
        std::cout << "\n";
    }

    /**
     * @brief pretty print queues
     * @note pass by value so the print uses a copy
     */
    template<typename T>
    static void printQueue(T q) 
    {
        while(!q.empty()) 
        {
            std::cout << q.top();
            if (q.size() != 1)
            {
                std::cout << ", ";
            } 
            q.pop();
        }
        std::cout << "\n";
    }
};

#endif // UTILS_H
