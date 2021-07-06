#include "Systematics.h"

#include <algorithm>

using namespace std;

Systematics::Systematics()
{
    m_systs = make_unique<vector<SystematicInfo*>>();
}

Systematics::~Systematics()
{
    for_each(content()->begin(), content()->end(), [](SystematicInfo* r) { delete r; });
    content()->clear();
}

void Systematics::add(const string &nm, const string &nmtex, eSystematicType tp, int col) const
{
    m_systs->emplace_back(new SystematicInfo(nm, nmtex, tp, col));
}
