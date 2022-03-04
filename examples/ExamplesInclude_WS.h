#include <string>
#include <vector>

using std::string;
using std::vector;

void test_ws_info(const string& filename);
int test_pulls(const string& filename, const string& outname);
int test_pulls_plot(const string& in, const string& out);
int test_ranking(const string& filename, const string& outname);
int test_ranking_plot(const string& in, const string& out);
void test_pseudodata(const vector<string>& path, const string& outfile);
void test_samplingdist(const string& filename, const std::string& tag);
void test_w2r(const string& filename, const string& postfit_result_file, const string& output);
void test_morphing(const string& filename);
