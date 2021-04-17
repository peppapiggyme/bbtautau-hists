#include "AutoBinningTool.h"

#include <fstream>
#include <sstream>

using namespace std;

#define BKGBKG (std::make_pair(eProcessType::BKG, eProcess::BKG))

// ----------------------------------------------------------------------------
// BinContent
// ----------------------------------------------------------------------------

std::size_t BinContent::size() const noexcept
{
    return m_size;
}

void BinContent::fill(const Config* c) noexcept
{
    vector<ProcessInfo*>* ps = c->processes->content();
    m_size = ps->front()->histogram->GetNbinsX();

    // total bkg
    data[BKGBKG] = std::make_pair(new double[m_size], new double[m_size]);
    for (size_t i = 0; i < m_size; ++i)
    {
        data[BKGBKG].first[i] = 0.;
        data[BKGBKG].second[i] = 0.;
    }
    for (const auto& p : *ps)
    {
        data[std::make_pair(p->type, p->process)] = std::move(getArrayFromContentAndError(p));
        if (p->type == eProcessType::BKG)
        {
            for (size_t i = 0; i < m_size; ++i)
            {
                data[BKGBKG].first[i] += data[std::make_pair(p->type, p->process)].first[i];
                data[BKGBKG].second[i] += std::pow(data[std::make_pair(p->type, p->process)].second[i], 2);
            }
        }
    }

    for (size_t i = 0; i < m_size; ++i)
    {
        data[BKGBKG].second[i] = std::sqrt(data[BKGBKG].second[i]);
    }
}

void BinContent::free() noexcept
{
    for (auto& d : data)
    {
        delete d.second.first;
        delete d.second.second;
    }
    data.clear();
}

std::pair<double*, double*> BinContent::getArrayFromContentAndError(ProcessInfo* p) const
{
    double* contents = new double[m_size];
    double* errors = new double[m_size];

    for (std::size_t i = 0; i < m_size; i++)
    {
        contents[i] = p->histogram->GetBinContent(i+1);
        errors[i] = p->histogram->GetBinError(i+1);
    }

    return std::move(std::make_pair(contents, errors));
}

void BinContent::showData() const
{
    for (auto& d : data)
    {
        Tools::println("%  %  [%(%), %(%), ..., %(%), %(%)]", 
            static_cast<int>(d.first.first), static_cast<int>(d.first.second), 
            d.second.first[0], d.second.second[0], 
            d.second.first[1], d.second.second[1],
            d.second.first[m_size-2], d.second.second[m_size-2],
            d.second.first[m_size-1], d.second.second[m_size-1]
        );
    }
}

bool BinContent::checkRangeXY(std::size_t x, std::size_t y) const
{
    if (x > y)
    {
        return false;
    }

    if (x > m_size - 1 || y > m_size - 1)
    {
        return false;
    }

    return true;
}

double BinContent::sumFromXtoY(double* arr, std::size_t x, std::size_t y) const
{
    if (!checkRangeXY(x, y)) 
    {
        return 0.;
    }

    if (x == y)
    {
        return arr[x];
    }

    double sum = 0.;

    for (size_t i = x; i <= y; ++i)
    {
        sum += arr[i];
    }

    return sum;
}

double BinContent::errSumFromXtoY(double* arr, std::size_t x, std::size_t y) const
{
    if (!checkRangeXY(x, y)) 
    {
        return 0.;
    }

    if (x == y)
    {
        return arr[x];
    }

    double sum2 = 0.;

    for (size_t i = x; i <= y; ++i)
    {
        sum2 += arr[i] * arr[i];
    }

    return std::sqrt(sum2);

}

double BinContent::significance(double fSig, double fBkg) const
{
    if (fSig < 0) fSig = 0;
    if (fBkg < 0) return 0;
    // Z = sqrt(2*(s+b)log(1+s/b)-2*s)
    return std::sqrt(2 * (fSig + fBkg) * std::log(1 + fSig / fBkg) - 2 * fSig);
} 

double BinContent::significanceBinByBin(std::priority_queue<std::size_t, std::vector<std::size_t>, std::greater<std::size_t>> pq) const
{
    auto itSig = std::find_if(data.begin(), data.end(), 
        [](const std::pair<std::pair<eProcessType, eProcess>, std::pair<double*, double*>>& d)
        { return d.first.first == eProcessType::SIG; });

    double fTotalSig2 = 0.;
    while (pq.size() > 1)
    {
        std::size_t p1 = pq.top(); pq.pop();
        std::size_t p2 = pq.top();
        double fSig = significance(sumFromXtoY(itSig->second.first, p1, p2), sumFromXtoY(data.at(BKGBKG).first, p1, p2));
        fTotalSig2 += fSig*fSig;
    }
    return std::sqrt(fTotalSig2);
}

bool BinContent::passCaseOneFromXtoY(double fBkg, double fBkgErr, double fReqFactor, std::size_t x, std::size_t y) const
{
    // if leading 2 main bkg pass this criteria, then it's okay
    // 2 -> total bkg and signal
    int nBkg = data.size() - 2;
    int nOkay = 2;
    for (const auto& d : data)
    {
        if (d.first.first == eProcessType::BKG)
        {
            double fErrOfThisProc = errSumFromXtoY(d.second.second, x, y) / sumFromXtoY(d.second.first, x, y);
            if (fErrOfThisProc > fReqFactor)
            {
                // Tools::println(" -> not good, Proc [%] Err [%] > %", static_cast<int>(d.first.second), fErrOfThisProc, fReqFactor);
                nBkg--;
            }
            // Tools::println(" <> good, Proc [%] Err [%] < %", static_cast<int>(d.first.second), fErrOfThisProc, fReqFactor);
        }
    }

    if (nBkg < nOkay)
    {
        return false;
    }

    double fTotalBkgErr = errSumFromXtoY(data.at(BKGBKG).second, x, y) / sumFromXtoY(data.at(BKGBKG).first, x, y);
    if (fTotalBkgErr > fBkgErr)
    {
        // Tools::println(" -> not good, Total Bkg Err [%] > %", fTotalBkgErr, fBkgErr);
        return false;
    }
    // Tools::println(" <> good, Total Bkg Err [%] < %", fTotalBkgErr, fBkgErr);

    double fTotalBkg = sumFromXtoY(data.at(BKGBKG).first, x, y);
    if (fTotalBkg < fBkg)
    {
        // Tools::println(" -> not good, Total Bkg [%] < %", fTotalBkg, fBkg);
        return false;
    }
    // Tools::println(" <> good, Total Bkg [%] > %", fTotalBkg, fBkg);
    return true;
}

bool BinContent::passCaseTwoFromXtoY(double fBkg, double fBkgErr, double fEffErr, double fReqFactor, std::size_t x, std::size_t y) const
{
    // if leading 2 main bkg pass this criteria, then it's okay
    // 2 -> total bkg and signal
    int nBkg = data.size() - 2;
    int nOkay = 2;

    auto itSig = std::find_if(data.begin(), data.end(), 
        [](const std::pair<std::pair<eProcessType, eProcess>, std::pair<double*, double*>>& d)
        { return d.first.first == eProcessType::SIG; });

    for (const auto& d : data)
    {
        if (d.first.first == eProcessType::BKG)
        {
            double fErrOfThisProc = errSumFromXtoY(d.second.second, x, y) / sumFromXtoY(d.second.first, x, y);
            if (fErrOfThisProc > fReqFactor)
            {
                // Tools::println(" -> not good, Proc [%] Err [%] > %", static_cast<int>(d.first.second), fErrOfThisProc, fReqFactor);
                nBkg--;
            }
            // Tools::println(" <> good, Proc [%] Err [%] < %", static_cast<int>(d.first.second), fErrOfThisProc, fReqFactor);
        }
    }

    if (nBkg < nOkay)
    {
        return false;
    }

    /// @note Hardcoded offset!
    constexpr double fOffset = 0.01;
    double fTotalSig = sumFromXtoY(itSig->second.first, x, y);
    double fTotalSigAll = sumFromXtoY(itSig->second.first, 0, size() - 1);
    double fTotalBkgErr = errSumFromXtoY(data.at(BKGBKG).second, x, y) / sumFromXtoY(data.at(BKGBKG).first, x, y);
    double fErr = fTotalSig > 0 ? std::min(fBkgErr, fOffset + fEffErr * (fTotalSig / fTotalSigAll)) : fBkgErr;

    if (fTotalBkgErr > fErr)
    {
        // Tools::println(" -> not good, Total Bkg Err [%] > %", fTotalBkgErr, fErr);
        return false;
    }
    // Tools::println(" <> good, Total Bkg Err [%] < %", fTotalBkgErr, fErr);

    double fTotalBkg = sumFromXtoY(data.at(BKGBKG).first, x, y);
    if (fTotalBkg < fBkg)
    {
        // Tools::println(" -> not good, Total Bkg [%] < %", fTotalBkg, fBkg);
        return false;
    }
    // Tools::println(" <> good, Total Bkg [%] > %", fTotalBkg, fBkg);
    return true;
}

// ----------------------------------------------------------------------------
// AutoBinningTool
// ----------------------------------------------------------------------------

AutoBinningTool::AutoBinningTool(const AutoBinningInfo* info)
    : m_info(info), m_binEdges(new std::vector<double>(m_info->n_bins + 1, 0.))
{
}

void AutoBinningTool::rebin(const Config* c) const 
{
    vector<double> y = *m_binEdges;
    vector<ProcessInfo*>* ps = c->processes->content();
    if (ps->front()->histogram->GetNbinsX() < (int)m_info->n_bins)
    {
        return;
    }

    for_each(ps->begin(), ps->end(), [this, &y, &c](ProcessInfo* p) {
        TH1* rebinned = p->histogram->Rebin(m_info->n_bins, p->histogram->GetName(), &y[0]);
        p->histogram = (TH1*)rebinned->Clone();
        for (auto& pp : p->systematic_histograms)
        {
            TH1* rebinned_pp = pp.second->Rebin(m_info->n_bins, p->histogram->GetName(), &y[0]);
            pp.second = (TH1*)rebinned_pp->Clone();
        } 
    });
}

// ----------------------------------------------------------------------------
// AutoBinningTool_v1
// ----------------------------------------------------------------------------

AutoBinningTool_v1::AutoBinningTool_v1(const AutoBinningInfo* info, BinningCriteria bc)
    : AutoBinningTool(info), m_bc(bc)
{
}

bool AutoBinningTool_v1::check(const Config* c) const
{
    if (!HistTool::check(c))
    {
        return false;
    }

    vector<ProcessInfo*>* ps = c->processes->content();

    int n_bkg_check = count_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return p->type == eProcessType::BKG; });

    int n_sig_check = count_if(ps->begin(), ps->end(), [](const ProcessInfo* p) {
        return p->type == eProcessType::SIG; });

    if (n_bkg_check < 1)
    {
        cerr << "FAIL: no background!\n";
        return false;
    }

    if (n_sig_check != 1) {
        cerr << "FAIL: must optimize with one signal!\n";
        return false;
    }

    return true;
}

void AutoBinningTool_v1::run(const Config* c) const
{
    m_binEdges->at(m_info->n_bins) = 1.;
    m_binEdges->at(0) = -1.;

    vector<std::size_t> curXs;

    BinContent* bc = new BinContent(c);

    bc->showData();

    bc->curX = bc->size() - 1;
    bc->curY = bc->size() - 1;

    // TODO: change this to binary searching
    while (bc->curX > 0)
    {
        bool bPass = false;

        // Choose binning criteria here
        switch (m_bc)
        {
        case BinningCriteria::CaseOne:
            bPass = bc->passCaseOneFromXtoY(m_info->min_bkg, m_info->min_mcstats, m_info->required_mcstats, bc->curX, bc->curY);
            break;
        case BinningCriteria::CaseTwo:
            bPass = bc->passCaseTwoFromXtoY(m_info->min_bkg, m_info->min_mcstats, m_info->eff_factor, m_info->required_mcstats, bc->curX, bc->curY);
            break;
        default: // e.g. None
            bPass = bc->passCaseTwoFromXtoY(m_info->min_bkg, m_info->min_mcstats, m_info->eff_factor, m_info->required_mcstats, bc->curX, bc->curY);
            break;
        }

        if (!bPass)
        {
            bc->curX--;
        }
        else
        {
            if (curXs.empty() || bc->curX != curXs.back())
            {
                // Tools::println("X: %  -  Y: %", bc->curX, bc->curY);
                curXs.push_back(bc->curX);
            } 
            else
            {
                bc->curX--;
            }
            // m_binEdges->at(i--) = -1. + 0.002 * (bc->curX + 1);
            bc->curY = bc->curX;
        }
    }

    std::priority_queue<std::size_t, std::vector<std::size_t>, std::greater<std::size_t>> pq;
    std::set<std::size_t> lookup;
    pq.push(0); pq.push(bc->size() - 1);

    auto pq_max(pq);
    for (auto x : curXs)
    {
        pq_max.push(x);
    }

    double sig_max = bc->significanceBinByBin(pq_max);
    Tools::println("Maximum one can get is %", sig_max);

    for (size_t i = 2; i <= m_info->n_bins; i++)
    {
        double bestSig = 0.;
        size_t bestX = 0;
        for (auto x : curXs) 
        {
            if (lookup.find(x) == lookup.end()) 
            {
                auto pq_copy(pq);
                pq_copy.push(x);
                double sig = bc->significanceBinByBin(pq_copy);
                if (sig > bestSig)
                {
                    bestSig = sig;
                    bestX = x;
                }
            }
        }
        // no repeated edge and ignore zero which has been pushed to pq already (see above)
        if (lookup.find(bestX) == lookup.end() && bestX != 0) {
            pq.push(bestX);
            lookup.insert(bestX);
        }
        if ((sig_max - bestSig) / sig_max < m_info->significance_delta)
        {
            break;
        }
    }

    double sig_final = bc->significanceBinByBin(pq);
    Tools::println("The optimal here is %", sig_final);
    Tools::println(" -> Difference is %", (sig_max - sig_final) / sig_max);

    // Tools::print_queue(pq);

    size_t i = 1;
    auto pq_forBinEdges = pq;
    pq_forBinEdges.pop(); // pop 0
    while(!pq_forBinEdges.empty() && i < m_info->n_bins)
    {
        m_binEdges->at(i++) = -1. + 0.002 * (pq_forBinEdges.top() + 1);
        pq_forBinEdges.pop();
    }

    vector<ProcessInfo*>* ps = c->processes->content();

    auto itSig = std::find_if(ps->begin(), ps->end(), 
        [](ProcessInfo* p) { return p->type == eProcessType::SIG; });

    ostringstream oss;
    oss << output_path << "/" 
        << "Binning_" 
        << Utils::signalTypeName((*itSig)->name) << ".txt";
    ofstream fout(oss.str());

    // Em.. WSMaker wants the opposite so
    std::priority_queue<size_t> pq_bin;
    auto pq_forBin = pq;
    while(!pq_forBin.empty())
    {
        int offset = 1;
        if (pq_forBin.size() == 1)
        {
            offset = 2; // WSMaker also wants overflow bin
        }
        pq_bin.push(pq_forBin.top() + offset); // bin number
        pq_forBin.pop();
    }

    Tools::print_queue(pq_bin);

    while(!pq_bin.empty())
    {
        fout << pq_bin.top();
        fout << " ";
        pq_bin.pop();
    }

    clog << "INFO: Binning saved in " << oss.str() << " for WSMaker" << '\n';
}

// ----------------------------------------------------------------------------
// AutoBinningTool_v2
// ----------------------------------------------------------------------------

AutoBinningTool_v2::AutoBinningTool_v2(const AutoBinningInfo* info, eProcess proc)
    : AutoBinningTool(info), m_proc(proc)
{
}

void AutoBinningTool_v2::run(const Config* c) const
{
    BinContent* bc = new BinContent(c);
    bc->showData();

    double* x = new double[m_info->n_bins];
    std::vector<double> y(m_info->n_bins + 1);

    for (size_t i = 0; i < m_info->n_bins; ++i)
    {
        x[i] = double(i+1) / (double)m_info->n_bins;
    }

    vector<ProcessInfo*>* ps = c->processes->content();
    auto it = std::find_if(ps->begin(), ps->end(), [this](ProcessInfo* p) { return p->process == m_proc; });
    TH1* h = (*it)->histogram;
    if (!h)
    {
        clog << "WARN: reference process not found, will use front\n";
        h = ps->front()->histogram;
    }
    else 
    {
        clog << "INFO: flattening based on " << static_cast<int>(m_proc) << '\n';
    }

    y[0] = h->GetXaxis()->GetBinLowEdge(1);
    h->GetQuantiles(m_info->n_bins, &y[1], x);
    *m_binEdges = std::move(y);
    delete x;
}
