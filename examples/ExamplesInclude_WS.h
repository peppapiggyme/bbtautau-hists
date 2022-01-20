#include <string>
#include <vector>

void test_ws_info(const std::string& filename);
int test_pulls(const std::string& filename, const std::string& outname);
int test_pulls_plot(const std::string& in, const std::string& out);
int test_ranking(const std::string& filename, const std::string& outname);
int test_ranking_plot(const std::string& in, const std::string& out);
void test_pseudodata(const std::vector<std::string>& path, const std::string& outfile);
void test_samplingdist(const std::string& filename);
void test_w2r(const std::string& filename, const std::string& postfit_result_file);