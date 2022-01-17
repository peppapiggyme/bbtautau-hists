

import ROOT as R

def qmu_muhat_lt_zero(mu, data, sig, bkg):
    result = 2.*mu*sum(sig) - 2.*sum([ni*R.TMath.Log(mu*si/bi+1.) for ni, si, bi in zip(data, sig, bkg)])
    return result if result > 0 else 0

# print(qmu_muhat_lt_zero(1., [1.], [1.], [1000.]))

# joint LTT and SLT
v9_data_last3 = [55.,5.,83.,10.,152.,20.]
v9_bkg_last3 = [59.738822,12.171817,104.014516,13.242181,151.758818,19.627705]
v9_sig_last3 = [8.027,2.614,7.534,2.725,7.356,2.611]

v13_data_last3 = [55.,5.,83.,10.,152.,20.]
v13_bkg_last3 = [58.024426,10.366922,101.863041,12.921958,149.54738,18.466209]
v13_sig_last3 = [7.932,2.585,7.446,2.696,7.271,2.583]

mu_0 = 0.4

print("v9", qmu_muhat_lt_zero(mu_0, v9_data_last3, v9_sig_last3, v9_bkg_last3))
print("v13", qmu_muhat_lt_zero(mu_0, v13_data_last3, v13_sig_last3, v13_bkg_last3))

print("v9", qmu_muhat_lt_zero(mu_0, [b for s, b in zip(v9_sig_last3, v9_bkg_last3)], v9_sig_last3, v9_bkg_last3))
print("v13", qmu_muhat_lt_zero(mu_0, [b for s, b in zip(v13_sig_last3, v13_bkg_last3)], v13_sig_last3, v13_bkg_last3))
