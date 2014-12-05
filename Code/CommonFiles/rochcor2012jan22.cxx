#include "rochcor2012jan22.h"
#include <TLorentzVector.h>


const float rochcor2012::netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};

const float rochcor2012::dcor_bf[8][8]={{-0.000019463,-0.000023303,0.000015430,0.000019712,0.000012505,0.000021642,0.000004119,-0.000076009},
					{0.000000128,0.000008017,0.000023615,0.000009331,0.000004833,0.000004227,0.000011144,-0.000035031},
					{-0.000028886,-0.000000872,0.000026070,-0.000004910,0.000005486,0.000006491,-0.000003761,-0.000013528},
					{-0.000059750,-0.000004919,0.000005412,-0.000003778,-0.000016824,-0.000013477,-0.000028910,-0.000083955},
					{-0.000118960,-0.000035253,-0.000000108,0.000007600,-0.000000344,-0.000001307,-0.000032807,-0.000076984},
					{-0.000030923,-0.000027221,0.000002016,0.000020695,0.000009981,0.000019933,-0.000005668,-0.000056440},
					{-0.000084168,-0.000033629,0.000009202,0.000009630,0.000009918,0.000027434,-0.000004246,-0.000050073},
					{-0.000065530,-0.000028725,0.000016278,0.000021019,0.000014359,0.000026534,0.000005743,0.000008129}};

const float rochcor2012::dcor_ma[8][8]={{0.000208866,0.000031300,-0.000020349,-0.000050529,-0.000053541,-0.000027905,0.000137902,0.000517059},
					{0.000495985,0.000113577,0.000007622,0.000000522,-0.000001945,0.000035598,0.000075060,0.000094852},
					{-0.000027405,-0.000011883,0.000029346,0.000047396,0.000051161,0.000017387,-0.000055606,-0.000242240},
					{-0.000403089,-0.000026173,0.000075550,0.000078408,0.000080089,0.000027961,-0.000070766,-0.000206802},
					{-0.000245396,0.000067964,0.000104941,0.000060063,0.000064881,0.000079890,0.000079182,-0.000031502},
					{-0.000070877,0.000004703,0.000016088,0.000012985,0.000024312,0.000049721,0.000105102,0.000108953},
					{0.000159998,-0.000005804,-0.000050732,-0.000032910,-0.000052791,-0.000051006,-0.000056044,0.000187873},
					{0.000038430,0.000015972,-0.000045405,-0.000058887,-0.000059324,-0.000113555,-0.000043767,0.000364713}};


const float rochcor2012::mcor_bf[8][8]={{-0.000071806,-0.000042738,-0.000005168,0.000008530,0.000006657,0.000001403,-0.000046476,-0.000077701},
					{-0.000079501,-0.000043963,0.000002809,0.000013418,0.000006740,-0.000009243,-0.000036724,-0.000053361},
					{-0.000086731,-0.000040736,-0.000007970,0.000002559,-0.000002886,-0.000003930,-0.000032835,-0.000108170},
					{-0.000078597,-0.000032426,-0.000000507,0.000008674,0.000011006,-0.000010700,-0.000049114,-0.000067586},
					{-0.000077383,-0.000039116,0.000001962,0.000007945,0.000006559,0.000002716,-0.000030194,-0.000072644},
					{-0.000047655,-0.000039136,0.000002394,0.000016047,0.000004600,-0.000002543,-0.000038468,-0.000069981},
					{-0.000057153,-0.000035963,-0.000001395,0.000013020,0.000013557,0.000000570,-0.000043211,-0.000073572},
					{-0.000078923,-0.000037795,-0.000002042,0.000014883,0.000010629,0.000002679,-0.000035537,-0.000058614}};
const float rochcor2012::mcor_ma[8][8]={{0.000282014,0.000062614,-0.000048363,-0.000058016,-0.000057209,-0.000077560,-0.000034640,0.000171967},
					{0.000146645,-0.000028641,-0.000046078,-0.000037986,-0.000040030,-0.000020762,0.000065515,0.000186970},
					{0.000013427,-0.000030694,0.000002951,0.000008824,0.000008922,0.000029666,0.000036620,-0.000070866},
					{0.000210925,0.000170459,0.000068899,0.000041447,0.000050304,0.000065185,0.000130217,0.000138370},
					{0.000216844,0.000137060,0.000075496,0.000056760,0.000056176,0.000063500,0.000115420,0.000237142},
					{-0.000361274,-0.000055447,0.000033291,0.000040406,0.000036132,0.000029008,0.000026157,0.000220532},
					{-0.000384866,-0.000080785,-0.000015095,-0.000009202,-0.000005052,-0.000003704,-0.000014772,-0.000031949},
					{0.000221685,0.000020988,-0.000065403,-0.000053498,-0.000045621,-0.000042911,-0.000154784,-0.000419276}};

//==========================================================================================================================
const float rochcor2012::dmavg[8][8]={{0.026037808,0.025249490,0.024935449,0.025273315,0.025256883,0.024945457,0.025334833,0.026298625},
				      {0.026134819,0.025244191,0.024973447,0.025380700,0.025404633,0.025003919,0.025175950,0.025845954},
				      {0.025777783,0.025134540,0.024971935,0.025387515,0.025403983,0.025023796,0.025113339,0.025517186},
				      {0.025631255,0.025206267,0.025040636,0.025380614,0.025432839,0.025050385,0.025179069,0.025668561},
				      {0.025804988,0.025336766,0.025030734,0.025293091,0.025362082,0.025067608,0.025326644,0.025864934},
				      {0.025737597,0.025235538,0.025007096,0.025349029,0.025375519,0.025014742,0.025273131,0.025892259},
				      {0.026018192,0.025225079,0.024968154,0.025297640,0.025267807,0.024928188,0.025138675,0.025986847},
				      {0.025919031,0.025259364,0.024935842,0.025269514,0.025302381,0.024888962,0.025176826,0.026070152}};
const float rochcor2012::dpavg[8][8]={{0.025794565,0.025216783,0.024989842,0.025384223,0.025380612,0.024993952,0.025132565,0.025608177},
				      {0.025425308,0.025022444,0.024939482,0.025396921,0.025419810,0.024972340,0.025105116,0.025768816},
				      {0.025922449,0.025250613,0.024938954,0.025365641,0.025362906,0.025003644,0.025286759,0.026077560},
				      {0.026311801,0.025282334,0.024897282,0.025294337,0.025330665,0.024980954,0.025352325,0.026137342},
				      {0.026267283,0.025269175,0.024905781,0.025272132,0.025255799,0.024937895,0.025208424,0.025977653},
				      {0.025948022,0.025231789,0.024975031,0.025351284,0.025359261,0.024896974,0.025130076,0.025782653},
				      {0.025886659,0.025302523,0.025028166,0.025368362,0.025365422,0.025031759,0.025249037,0.025750145},
				      {0.025887077,0.025243519,0.025030559,0.025399390,0.025416798,0.025092332,0.025264352,0.025626461}};

const float rochcor2012::mmavg[8][8]={{0.026058156,0.025288242,0.024914332,0.025234715,0.025215595,0.024929635,0.025270115,0.025953079},
				      {0.025884172,0.025161363,0.024956277,0.025306438,0.025303837,0.024965949,0.025254987,0.025902839},
				      {0.025917520,0.025163348,0.024993308,0.025319856,0.025347663,0.025035577,0.025274874,0.025796417},
				      {0.026117465,0.025390911,0.025062392,0.025311546,0.025316478,0.025039940,0.025369449,0.025927631},
				      {0.026077926,0.025368480,0.025004913,0.025266120,0.025282883,0.025055789,0.025296928,0.026123998},
				      {0.025568735,0.025178538,0.025031160,0.025364352,0.025371157,0.025003405,0.025254490,0.025926064},
				      {0.025522195,0.025161997,0.024978428,0.025258552,0.025309075,0.024950539,0.025195544,0.025826596},
				      {0.026082824,0.025253782,0.024939655,0.025269799,0.025249572,0.024941604,0.025109541,0.025501534}};
const float rochcor2012::mpavg[8][8]={{0.025715677,0.025167430,0.025012492,0.025360102,0.025379251,0.025026302,0.025263846,0.025867599},
				      {0.025773304,0.025212873,0.024995786,0.025398794,0.025404635,0.025042453,0.025220045,0.025789688},
				      {0.025956906,0.025259180,0.024997460,0.025345088,0.025376761,0.024991251,0.025196061,0.025988617},
				      {0.025820694,0.025128369,0.024896937,0.025237798,0.025295188,0.024883958,0.025165118,0.025862797},
				      {0.025796839,0.025146159,0.024926421,0.025250195,0.025246977,0.024944807,0.025171534,0.025776059},
				      {0.026171851,0.025266560,0.024932881,0.025320040,0.025318938,0.024982888,0.025231810,0.025805723},
				      {0.026380484,0.025369597,0.025003962,0.025313980,0.025316192,0.025017893,0.025234826,0.026008156},
				      {0.025769622,0.025252647,0.025040967,0.025371246,0.025346839,0.025008565,0.025377072,0.026334648}};
//==========================================================================================================================
const float rochcor2012::dcor_bfer[8][8]={{0.000043517,0.000022850,0.000018965,0.000017370,0.000017296,0.000019165,0.000022475,0.000041328},
					  {0.000039740,0.000021864,0.000019136,0.000017649,0.000017492,0.000019318,0.000021809,0.000039973},
					  {0.000040881,0.000022519,0.000019012,0.000017486,0.000017326,0.000018985,0.000021917,0.000039925},
					  {0.000039796,0.000022106,0.000018885,0.000017362,0.000017267,0.000018970,0.000021768,0.000039489},
					  {0.000041281,0.000022815,0.000019045,0.000017438,0.000017376,0.000019031,0.000022383,0.000042126},
					  {0.000040000,0.000021836,0.000018868,0.000017763,0.000017724,0.000019005,0.000022178,0.000039864},
					  {0.000040347,0.000022329,0.000018999,0.000017338,0.000017632,0.000018900,0.000021596,0.000039967},
					  {0.000040626,0.000022327,0.000018884,0.000017354,0.000017223,0.000019075,0.000021874,0.000039248}};
const float rochcor2012::dcor_maer[8][8]={{0.000043517,0.000022850,0.000018965,0.000017370,0.000017296,0.000019165,0.000022475,0.000041328},
					  {0.000039740,0.000021864,0.000019136,0.000017649,0.000017492,0.000019318,0.000021809,0.000039973},
					  {0.000040881,0.000022519,0.000019012,0.000017486,0.000017326,0.000018985,0.000021917,0.000039925},
					  {0.000039796,0.000022106,0.000018885,0.000017362,0.000017267,0.000018970,0.000021768,0.000039489},
					  {0.000041281,0.000022815,0.000019045,0.000017438,0.000017376,0.000019031,0.000022383,0.000042126},
					  {0.000040000,0.000021836,0.000018868,0.000017763,0.000017724,0.000019005,0.000022178,0.000039864},
					  {0.000040347,0.000022329,0.000018999,0.000017338,0.000017632,0.000018900,0.000021596,0.000039967},
					  {0.000040626,0.000022327,0.000018884,0.000017354,0.000017223,0.000019075,0.000021874,0.000039248}};

const float rochcor2012::mcor_bfer[8][8]={{0.000047815,0.000025172,0.000021069,0.000019264,0.000019165,0.000021343,0.000024927,0.000045633},
					  {0.000044371,0.000024478,0.000021293,0.000019543,0.000019431,0.000021466,0.000024209,0.000044329},
					  {0.000045310,0.000025143,0.000021216,0.000019483,0.000019279,0.000021192,0.000024523,0.000044429},
					  {0.000044504,0.000024779,0.000021058,0.000019274,0.000019184,0.000021143,0.000024413,0.000044281},
					  {0.000045494,0.000025065,0.000021169,0.000019306,0.000019304,0.000021131,0.000024836,0.000046798},
					  {0.000044599,0.000024474,0.000020983,0.000019729,0.000019618,0.000021050,0.000024673,0.000044435},
					  {0.000045441,0.000024933,0.000021162,0.000019292,0.000019560,0.000021106,0.000024186,0.000044568},
					  {0.000045326,0.000024917,0.000021085,0.000019279,0.000019166,0.000021262,0.000024469,0.000043777}};
const float rochcor2012::mcor_maer[8][8]={{0.000047815,0.000025172,0.000021069,0.000019264,0.000019165,0.000021343,0.000024927,0.000045633},
					  {0.000044371,0.000024478,0.000021293,0.000019543,0.000019431,0.000021466,0.000024209,0.000044329},
					  {0.000045310,0.000025143,0.000021216,0.000019483,0.000019279,0.000021192,0.000024523,0.000044429},
					  {0.000044504,0.000024779,0.000021058,0.000019274,0.000019184,0.000021143,0.000024413,0.000044281},
					  {0.000045494,0.000025065,0.000021169,0.000019306,0.000019304,0.000021131,0.000024836,0.000046798},
					  {0.000044599,0.000024474,0.000020983,0.000019729,0.000019618,0.000021050,0.000024673,0.000044435},
					  {0.000045441,0.000024933,0.000021162,0.000019292,0.000019560,0.000021106,0.000024186,0.000044568},
					  {0.000045326,0.000024917,0.000021085,0.000019279,0.000019166,0.000021262,0.000024469,0.000043777}};

const float rochcor2012::sf[8] = {0.0153729,0.0103115,0.00701322,0.00472529,0.00460413,0.00700919,0.00967325,0.0147967};
const float rochcor2012:: sfer[8] = {0.00026422,0.000131926,0.000158204,0.000141117,0.000152415,0.000145095,0.000128231,0.000244732};

const float rochcor2012::gsf[8] = {0.998806,0.999243,0.999688,1.00009,1.00016,0.999799,0.999237,0.998916};
const float rochcor2012::gsfer[8] = {0.000170426,6.26483e-05,4.65517e-05,3.69658e-05,3.70135e-05,4.67504e-05,6.13789e-05,0.000168485};

//===============================================================================================
//parameters for Z pt correction

const float rochcor2012::ptlow[85] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5,
				      6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
				      10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
				      15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5,
				      20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5,
				      25.0, 26.0, 27.0, 28.0, 29.0,
				      30.0, 32.0, 34.0, 36.0, 38.0,
				      40.0, 44.0, 48.0, 52.0, 56.0,
				      60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0,
				      100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 175.0,
				      200.0, 250.0, 350.0, 500.0, 1000.0};

//int nptbins( sizeof(ptlow)/sizeof(float) - 1 );

const float rochcor2012::zptscl[84] = {1.49177,1.45654,1.36283,1.28569,1.2418,1.12336,1.10416,1.08731,0.994051,0.96532,
				       0.94427,0.932725,0.918082,0.899665,0.898398,0.927687,0.908047,0.892392,0.924027,0.945895,
				       0.937149,0.923983,0.923387,0.955362,0.947812,0.962943,0.948781,0.961555,0.95222,0.999207,
				       0.973884,0.993013,0.953487,0.951402,0.985583,0.986603,0.981388,1.00022,1.0294,0.964748,
				       0.974592,1.01546,0.992343,1.00101,0.990866,0.98982,1.02924,1.02265,0.967695,1.02411,
				       0.97331,1.01052,1.01561,0.992594,0.976504,1.01205,0.981111,1.00078,1.02078,1.00719,
				       1.0099,1.02865,1.03845,1.03254,1.09815,1.10263,1.06302,1.0725,1.14703,1.10574,
				       1.13911,1.16947,1.1709,1.11413,1.28793,1.18953,1.20212,1.18112,1.25471,1.15329,
				       1.14276,1.17223,1.09173,2.00229};

const float rochcor2012::zptscler[84] = {0.0270027,0.0154334,0.0115338,0.00958085,0.0084683,0.00736665,0.0069567,0.00671434,
					 0.00617693,0.00601943,0.00594735,0.00594569,0.00594903,0.00595495,0.00608115,0.00633704,
					 0.0063916,0.0064468,0.00678106,0.00706769,0.00717517,0.00727958,0.00747182,0.00785544,
					 0.00798754,0.00828787,0.00839147,0.00865826,0.00876775,0.00933276,0.00935768,0.0097289,
					 0.00962058,0.00983828,0.0103044,0.0104871,0.0106575,0.0110388,0.0114986,0.0111494,
					 0.0115202,0.0121059,0.0121345,0.0124923,0.0125972,0.0128401,0.0134519,0.0136279,
					 0.0133414,0.014186,0.00992195,0.0105984,0.0109484,0.0111756,0.0114579,0.00870013,
					 0.00904749,0.00970734,0.0104583,0.0109818,0.00837852,0.00939894,0.010415,0.0113433,
					 0.013007,0.0128788,0.0140174,0.0156993,0.0181717,0.019765,0.0222326,0.0249408,
					 0.0272806,0.0211706,0.0278087,0.0306654,0.0361387,0.041327,0.0341513,0.0440116,
					 0.0473006,0.0680212,0.149162,0.56279};

rochcor2012::~rochcor2012(){
}

rochcor2012::rochcor2012(){
  
  eran.SetSeed(123456);
  sran.SetSeed(1234);
  
  gscler_mc_dev=0;
  gscler_da_dev=0;

  for(int i=0; i<8; ++i){
      for(int j=0; j<8; ++j){
          mptsys_mc_dm[i][j]=0;
          mptsys_mc_da[i][j]=0;
          mptsys_da_dm[i][j]=0;
          mptsys_da_da[i][j]=0;
      }
  }

}

rochcor2012::rochcor2012(int seed){
  eran.SetSeed(123456);
  sran.SetSeed(seed);

  gscler_mc_dev=sran.Gaus(0.0, 1.0);
  gscler_da_dev=sran.Gaus(0.0, 1.0);

  for(int i=0; i<8; ++i){
      for(int j=0; j<8; ++j){
          mptsys_mc_dm[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_mc_da[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_da_dm[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_da_da[i][j]=sran.Gaus(0.0, 1.0);
      }
  }
}

void rochcor2012::momcor_mc( TLorentzVector& mu, float charge, int runopt, float& qter){
  
  //sysdev == num : deviation = num

  float ptmu = mu.Pt();
  float muphi = mu.Phi();
  float mueta = mu.Eta(); // same with mu.Eta() in Root

  float px = mu.Px();
  float py = mu.Py();
  float pz = mu.Pz();
  float e = mu.E();

  int mu_phibin = phibin(muphi);
  int mu_etabin = etabin(mueta);
  
  if(mu_phibin>=0 && mu_etabin>=0){
    
    float Mf = (mcor_bf[mu_phibin][mu_etabin] + mptsys_mc_dm[mu_phibin][mu_etabin]*mcor_bfer[mu_phibin][mu_etabin])/(mpavg[mu_phibin][mu_etabin]+mmavg[mu_phibin][mu_etabin]);
    float Af = ((mcor_ma[mu_phibin][mu_etabin]+mptsys_mc_da[mu_phibin][mu_etabin]*mcor_maer[mu_phibin][mu_etabin]) - Mf*(mpavg[mu_phibin][mu_etabin]-mmavg[mu_phibin][mu_etabin]));     
    
    float cor = 1.0/(1.0 + 2.0*Mf + charge*Af*ptmu);
    
    //for the momentum tuning - eta,phi,Q correction
    px *= cor;
    py *= cor;
    pz *= cor;
    e  *= cor;
    
    float gscler = mgscl_stat;    
    float gscl = (genm_smr/mrecm);
    
    px *= (gscl + gscler_mc_dev*gscler);
    py *= (gscl + gscler_mc_dev*gscler);
    pz *= (gscl + gscler_mc_dev*gscler);
    e  *= (gscl + gscler_mc_dev*gscler);
    
    float momscl = sqrt(px*px + py*py)/ptmu;
    
    float tune = gsf[mu_etabin]*eran.Gaus(1.0,sf[mu_etabin]);
    
    px *= (tune); 
    py *= (tune);  
    pz *= (tune);  
    e  *= (tune);   
    
    qter *= sqrt(momscl*momscl + (1.0-tune)*(1.0-tune));
  }
  
  mu.SetPxPyPzE(px,py,pz,e);
  
}


void rochcor2012::momcor_data( TLorentzVector& mu, float charge, int runopt, float& qter){
  
  float ptmu = mu.Pt();

  float muphi = mu.Phi();
  float mueta = mu.Eta(); // same with mu.Eta() in Root

  float px = mu.Px();
  float py = mu.Py();
  float pz = mu.Pz();
  float e = mu.E();
  
  int mu_phibin = phibin(muphi);
  int mu_etabin = etabin(mueta);

  if(mu_phibin>=0 && mu_etabin>=0){

    float Mf = (dcor_bf[mu_phibin][mu_etabin]+mptsys_da_dm[mu_phibin][mu_etabin]*dcor_bfer[mu_phibin][mu_etabin])/(dpavg[mu_phibin][mu_etabin]+dmavg[mu_phibin][mu_etabin]);
    float Af = ((dcor_ma[mu_phibin][mu_etabin]+mptsys_da_da[mu_phibin][mu_etabin]*dcor_maer[mu_phibin][mu_etabin]) - Mf*(dpavg[mu_phibin][mu_etabin]-dmavg[mu_phibin][mu_etabin]));     
    
    float cor = 1.0/(1.0 + 2.0*Mf + charge*Af*ptmu);
    
    px *= cor;
    py *= cor;
    pz *= cor;
    e  *= cor;
    
    //after Z pt correction
    float gscler = dgscl_stat;
    float gscl = (genm_smr/drecm);
    
    px *= (gscl + gscler_da_dev*gscler);
    py *= (gscl + gscler_da_dev*gscler);
    pz *= (gscl + gscler_da_dev*gscler);
    e  *= (gscl + gscler_da_dev*gscler);
    
    float momscl = sqrt(px*px + py*py)/ptmu;
    qter *= momscl;
    
  }
  
  mu.SetPxPyPzE(px,py,pz,e);
  
}

int rochcor2012::phibin(float phi){
  
  int nphibin = -1;
  
  for(int i=0; i<8; i++){
    if(-1*pi+(2.0*pi/8.0)*i <= phi && -1*pi+(2.0*pi/8.0)*(i+1) > phi){
      nphibin = i;
      break;
    }
  }
  
  return nphibin;
}

int rochcor2012::etabin(float eta){

  int nbin = -1;
  
  for(int i=0; i<8; i++){
    if(netabin[i] <= eta && netabin[i+1] > eta){
      nbin = i;
      break;
    }
  }
  
  return nbin;
}

float rochcor2012::zptcor(float gzpt) {
  int ibin( 0 );
  
  // mcptscl[] = 84 bins: [0] and [83] are the underflow and overflow
  if ( gzpt > ptlow[nptbins] ) return nptbins-1;
  if ( gzpt < ptlow[0      ] ) return 0;
  
  for ( int i=0; i<nptbins; ++i ) {
    if ( gzpt>=ptlow[i] && gzpt<ptlow[i+1] ) { ibin=i; break; }
  }

  float zptwt = zptscl[ibin];

  return zptwt;
}