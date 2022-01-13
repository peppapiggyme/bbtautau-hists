#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

#include <memory>
#include <vector>
#include <string>

using std::string;
using std::vector;
using std::unique_ptr;

// EUMS
enum class eSystematicType { TwoSide, OneSide };


// ELEM
class SystematicInfo
{
public:
    SystematicInfo(const string& nm, const string& nmtex, eSystematicType tp, int col=2, bool smth=false) noexcept
        : name(nm), name_tex(nmtex), type(tp), color(col), smooth(smth) {}

public:
    string name;
    string name_tex;
    eSystematicType type;
    int color;
    bool smooth;  // only monotonic smoothing, as in WSMaker
};


// COLLECTION
class Systematics
{
public:
    Systematics();
    ~Systematics();
    Systematics(Systematics& rs) = delete;
    Systematics& operator=(Systematics& rs) = delete;

public:
    void add(const string& nm, const string& nmtex, eSystematicType tp, int col=2, bool smth=false) const;
    inline vector<SystematicInfo*>* content() const { return m_systs.get(); }

private:
    unique_ptr<vector<SystematicInfo*>> m_systs;
};

#endif // SYSTEMATICS_H