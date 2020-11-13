//=====================================================================================
//---------------------------Declaration by the Developer------------------------------
//=====================================================================================
/*
	This program(library) is wrote for scientific reasearch usage by Yiping Hao who began his graduate career since 2019.
	When the code was writen, the author was working for Dalian Institute of Chemical Physics in researcher Donghui Zhang's group.
	Its usages are beyond the research fields of Donghui Zhang's group. And it's not the crucial codes of Donghui Zhang's group.
    Thus the author gives up its copyright.
    Anyone who uses it in his/her research should express gratitude in his/her corresponding published papers.
    The author is not responsible for any adverse consequences caused by the use of this code. 
	The first release was on 2020.10.30.
	The author's email addresses are Hyper@mail.ustc.edu.cn and Hyper@dicp.ac.cn. His QQ number is 649405039.
	All the communications refer to this program(library) will be welcomed.
	If a viewer wants to modifies this code, there is no necessity informing Yiping Hao.
	However, she or he must rename the program&code and refresh the Declaration with its own name.
	Yiping Hao should not be blamed for problems caused by the new version.
	*********************To your valour, my sword, and our victory together. Long may the sun shine!**********************
*/
//=====================================================================================
//--------------------------------About the Library------------------------------------
//=====================================================================================
/*
	This Library is about processes that have general usage in computional chemsitry.
    Funcitions in this library almost are lightweight and never with complex mathematical computation.
    CONTENTS
	Part One: Database of Elements and Isotopes' Masses.
*/
struct mass_list
{
    int proton;
    int neutron;
    double mass_of_atom;//relative atomic mass unit
    double mass_of_nucleus;//atomic unit, mass of atom - mass of static eletrons (without reckoning binding energy)
    double RAM_of_element;//keep the same for the same element//relative atomic mass
    double abundency;
    char symbol[4];
    char name[20];
};
#define _RAM_LIST_LENGTH_ 322;
static const double RAM_tu_AU = 1822.88761;
static struct mass_list list[_RAM_LIST_LENGTH_];
static int list_if_inilialized = 0;
//=====================================================================================
//-----------------------------------Declarations--------------------------------------
//=====================================================================================
//-----------------Part One: Database of Elements and Isotopes' Masses-----------------
static void value_mass_list(struct mass_list *list);
double AU_UNIT_MASS(int amount_proton, int amount_neutron);
static void STR_NAME_COPY(char*to, char*from);
//=====================================================================================
//--------------------------------------Library----------------------------------------
//=====================================================================================
static void value_mass_list(struct mass_list *list)
{
    /*
    This table lists the mass and percent natural abundance for the stable nuclides.  
    The mass of the longest lived isotope is given for elements without a stable nuclide.    
    The isotopic mass data is from G. Audi, A. H. Wapstra Nucl. Phys A. 1993, 565, 1-65 and G. Audi, A. H. Wapstra Nucl. Phys A. 1995, 595, 409-480.  
    The percent natural abundance data is from the 1997 report of the IUPAC Subcommittee for Isotopic Abundance Measurements by K.J.R. Rosman, P.D.P. Taylor Pure Appl. Chem. 1999, 71, 1593-1607. 
    */
    //==================================neutron==================================
    (list[0]).proton = 0;
    (list[0]).neutron = 1;
    (list[0]).mass = 1.0086654;
    (list[0]).mass_of_nucleus = RAM_tu_AU * 1.0086654;
    (list[0]).RAM_of_element = 1.0086654;
    (list[0]).abundency = 1.0000000;
    STR_NAME_COPY((list[0]).symbol, "n");
    STR_NAME_COPY((list[0]).name, "Neutron");
    //===========================================================================
    (list[1]).proton = 1;
    (list[1]).neutron = 0;
    (list[1]).mass_of_atom = 1.007825;
    (list[1]).mass_of_nucleus = 1836.15170554825;
    (list[1]).RAM_of_element = 1.0079;
    (list[1]).abundency = 0.999885;
    STR_NAME_COPY((list[1]).symbol, "H");
    STR_NAME_COPY((list[1]).name, "Hydrogen");
    //===========================================================================
    (list[2]).proton = 1;
    (list[2]).neutron = 1;
    (list[2]).mass_of_atom = 2.014102;
    (list[2]).mass_of_nucleus = 3670.4815810762198;
    (list[2]).RAM_of_element = 1.0079;
    (list[2]).abundency = 0.000115;
    STR_NAME_COPY((list[2]).symbol, "H");
    STR_NAME_COPY((list[2]).name, "Deuterium");
    //===========================================================================
    (list[3]).proton = 1;
    (list[3]).neutron = 2;
    (list[3]).mass_of_atom = 3.016049;
    (list[3]).mass_of_nucleus = 5496.9183532528905;
    (list[3]).RAM_of_element = 1.0079;
    (list[3]).abundency = 0.0;
    STR_NAME_COPY((list[3]).symbol, "H");
    STR_NAME_COPY((list[3]).name, "Tritium");
    //===========================================================================
    (list[4]).proton = 2;
    (list[4]).neutron = 1;
    (list[4]).mass_of_atom = 3.016029;
    (list[4]).mass_of_nucleus = 5495.88189550069;
    (list[4]).RAM_of_element = 4.0026;
    (list[4]).abundency = 1.37e-06;
    STR_NAME_COPY((list[4]).symbol, "He");
    STR_NAME_COPY((list[4]).name, "Helium");
    //===========================================================================
    (list[5]).proton = 2;
    (list[5]).neutron = 2;
    (list[5]).mass_of_atom = 4.002603;
    (list[5]).mass_of_nucleus = 7294.295416448829;
    (list[5]).RAM_of_element = 4.0026;
    (list[5]).abundency = 0.99999863;
    STR_NAME_COPY((list[5]).symbol, "He");
    STR_NAME_COPY((list[5]).name, "Helium");
    //===========================================================================
    (list[6]).proton = 3;
    (list[6]).neutron = 3;
    (list[6]).mass_of_atom = 6.015122;
    (list[6]).mass_of_nucleus = 10961.89136643842;
    (list[6]).RAM_of_element = 6.941;
    (list[6]).abundency = 0.0759;
    STR_NAME_COPY((list[6]).symbol, "Li");
    STR_NAME_COPY((list[6]).name, "Lithium");
    //===========================================================================
    (list[7]).proton = 3;
    (list[7]).neutron = 4;
    (list[7]).mass_of_atom = 7.016004;
    (list[7]).mass_of_nucleus = 12786.38676331044;
    (list[7]).RAM_of_element = 6.941;
    (list[7]).abundency = 0.9240999999999999;
    STR_NAME_COPY((list[7]).symbol, "Li");
    STR_NAME_COPY((list[7]).name, "Lithium");
    //===========================================================================
    (list[8]).proton = 4;
    (list[8]).neutron = 5;
    (list[8]).mass_of_atom = 9.012182;
    (list[8]).mass_of_nucleus = 16424.194906865017;
    (list[8]).RAM_of_element = 9.0122;
    (list[8]).abundency = 1.0;
    STR_NAME_COPY((list[8]).symbol, "Be");
    STR_NAME_COPY((list[8]).name, "Beryllium");
    //===========================================================================
    (list[9]).proton = 5;
    (list[9]).neutron = 5;
    (list[9]).mass_of_atom = 10.012937;
    (list[9]).mass_of_nucleus = 18247.45879701057;
    (list[9]).RAM_of_element = 10.811;
    (list[9]).abundency = 0.19899999999999998;
    STR_NAME_COPY((list[9]).symbol, "B");
    STR_NAME_COPY((list[9]).name, "Boron");
    //===========================================================================
    (list[10]).proton = 5;
    (list[10]).neutron = 6;
    (list[10]).mass_of_atom = 11.009305;
    (list[10]).mass_of_nucleus = 20063.72567921105;
    (list[10]).RAM_of_element = 10.811;
    (list[10]).abundency = 0.8009999999999999;
    STR_NAME_COPY((list[10]).symbol, "B");
    STR_NAME_COPY((list[10]).name, "Boron");
    //===========================================================================
    (list[11]).proton = 6;
    (list[11]).neutron = 6;
    (list[11]).mass_of_atom = 12.000000;
    (list[11]).mass_of_nucleus = 21868.65132;
    (list[11]).RAM_of_element = 12.0107;
    (list[11]).abundency = 0.9893000000000001;
    STR_NAME_COPY((list[11]).symbol, "C");
    STR_NAME_COPY((list[11]).name, "Carbon");
    //===========================================================================
    (list[12]).proton = 6;
    (list[12]).neutron = 7;
    (list[12]).mass_of_atom = 13.003355;
    (list[12]).mass_of_nucleus = 23697.654717931553;
    (list[12]).RAM_of_element = 12.0107;
    (list[12]).abundency = 0.010700000000000001;
    STR_NAME_COPY((list[12]).symbol, "C");
    STR_NAME_COPY((list[12]).name, "Carbon");
    //===========================================================================
    (list[13]).proton = 6;
    (list[13]).neutron = 8;
    (list[13]).mass_of_atom = 14.003242;
    (list[13]).mass_of_nucleus = 25520.33634163162;
    (list[13]).RAM_of_element = 12.0107;
    (list[13]).abundency = 0.0;
    STR_NAME_COPY((list[13]).symbol, "C");
    STR_NAME_COPY((list[13]).name, "Carbon");
    //===========================================================================
    (list[14]).proton = 7;
    (list[14]).neutron = 7;
    (list[14]).mass_of_atom = 14.003074;
    (list[14]).mass_of_nucleus = 25519.03009651314;
    (list[14]).RAM_of_element = 14.0067;
    (list[14]).abundency = 0.9963200000000001;
    STR_NAME_COPY((list[14]).symbol, "N");
    STR_NAME_COPY((list[14]).name, "Nitrogen");
    //===========================================================================
    (list[15]).proton = 7;
    (list[15]).neutron = 8;
    (list[15]).mass_of_atom = 15.000109;
    (list[15]).mass_of_nucleus = 27336.51284474949;
    (list[15]).RAM_of_element = 14.0067;
    (list[15]).abundency = 0.00368;
    STR_NAME_COPY((list[15]).symbol, "N");
    STR_NAME_COPY((list[15]).name, "Nitrogen");
    //===========================================================================
    (list[16]).proton = 8;
    (list[16]).neutron = 8;
    (list[16]).mass_of_atom = 15.994915;
    (list[16]).mass_of_nucleus = 29148.93237650315;
    (list[16]).RAM_of_element = 15.9994;
    (list[16]).abundency = 0.9975700000000001;
    STR_NAME_COPY((list[16]).symbol, "O");
    STR_NAME_COPY((list[16]).name, "Oxygen");
    //===========================================================================
    (list[17]).proton = 8;
    (list[17]).neutron = 9;
    (list[17]).mass_of_atom = 16.999132;
    (list[17]).mass_of_nucleus = 30979.50710355452;
    (list[17]).RAM_of_element = 15.9994;
    (list[17]).abundency = 0.00037999999999999997;
    STR_NAME_COPY((list[17]).symbol, "O");
    STR_NAME_COPY((list[17]).name, "Oxygen");
    //===========================================================================
    (list[18]).proton = 8;
    (list[18]).neutron = 10;
    (list[18]).mass_of_atom = 17.999160;
    (list[18]).mass_of_nucleus = 32802.4457544076;
    (list[18]).RAM_of_element = 15.9994;
    (list[18]).abundency = 0.0020499999999999997;
    STR_NAME_COPY((list[18]).symbol, "O");
    STR_NAME_COPY((list[18]).name, "Oxygen");
    //===========================================================================
    (list[19]).proton = 9;
    (list[19]).neutron = 10;
    (list[19]).mass_of_atom = 18.998403;
    (list[19]).mass_of_nucleus = 34622.95343848683;
    (list[19]).RAM_of_element = 18.9984;
    (list[19]).abundency = 1.0;
    STR_NAME_COPY((list[19]).symbol, "F");
    STR_NAME_COPY((list[19]).name, "Fluorine");
    //===========================================================================
    (list[20]).proton = 10;
    (list[20]).neutron = 10;
    (list[20]).mass_of_atom = 19.992440;
    (list[20]).mass_of_nucleus = 36433.971169668395;
    (list[20]).RAM_of_element = 20.1797;
    (list[20]).abundency = 0.9048;
    STR_NAME_COPY((list[20]).symbol, "Ne");
    STR_NAME_COPY((list[20]).name, "Neon");
    //===========================================================================
    (list[21]).proton = 10;
    (list[21]).neutron = 11;
    (list[21]).mass_of_atom = 20.993847;
    (list[21]).mass_of_nucleus = 38259.423582535666;
    (list[21]).RAM_of_element = 20.1797;
    (list[21]).abundency = 0.0027;
    STR_NAME_COPY((list[21]).symbol, "Ne");
    STR_NAME_COPY((list[21]).name, "Neon");
    //===========================================================================
    (list[22]).proton = 10;
    (list[22]).neutron = 12;
    (list[22]).mass_of_atom = 21.991386;
    (list[22]).mass_of_nucleus = 40077.82506612746;
    (list[22]).RAM_of_element = 20.1797;
    (list[22]).abundency = 0.0925;
    STR_NAME_COPY((list[22]).symbol, "Ne");
    STR_NAME_COPY((list[22]).name, "Neon");
    //===========================================================================
    (list[23]).proton = 11;
    (list[23]).neutron = 12;
    (list[23]).mass_of_atom = 22.989770;
    (list[23]).mass_of_nucleus = 41896.7668897497;
    (list[23]).RAM_of_element = 22.9897;
    (list[23]).abundency = 1.0;
    STR_NAME_COPY((list[23]).symbol, "Na");
    STR_NAME_COPY((list[23]).name, "Sodium");
    //===========================================================================
    (list[24]).proton = 12;
    (list[24]).neutron = 12;
    (list[24]).mass_of_atom = 23.985042;
    (list[24]).mass_of_nucleus = 43710.03588712962;
    (list[24]).RAM_of_element = 24.305;
    (list[24]).abundency = 0.7898999999999999;
    STR_NAME_COPY((list[24]).symbol, "Mg");
    STR_NAME_COPY((list[24]).name, "Magnesium");
    //===========================================================================
    (list[25]).proton = 12;
    (list[25]).neutron = 13;
    (list[25]).mass_of_atom = 24.985837;
    (list[25]).mass_of_nucleus = 45534.37269277957;
    (list[25]).RAM_of_element = 24.305;
    (list[25]).abundency = 0.1;
    STR_NAME_COPY((list[25]).symbol, "Mg");
    STR_NAME_COPY((list[25]).name, "Magnesium");
    //===========================================================================
    (list[26]).proton = 12;
    (list[26]).neutron = 14;
    (list[26]).mass_of_atom = 25.982593;
    (list[26]).mass_of_nucleus = 47351.34685537273;
    (list[26]).RAM_of_element = 24.305;
    (list[26]).abundency = 0.1101;
    STR_NAME_COPY((list[26]).symbol, "Mg");
    STR_NAME_COPY((list[26]).name, "Magnesium");
    //===========================================================================
    (list[27]).proton = 13;
    (list[27]).neutron = 14;
    (list[27]).mass_of_atom = 26.981538;
    (list[27]).mass_of_nucleus = 49171.31131894418;
    (list[27]).RAM_of_element = 26.9815;
    (list[27]).abundency = 1.0;
    STR_NAME_COPY((list[27]).symbol, "Al");
    STR_NAME_COPY((list[27]).name, "Aluminum");
    //===========================================================================
    (list[28]).proton = 14;
    (list[28]).neutron = 14;
    (list[28]).mass_of_atom = 27.976927;
    (list[28]).mass_of_nucleus = 50984.79359417447;
    (list[28]).RAM_of_element = 28.0855;
    (list[28]).abundency = 0.9222969999999999;
    STR_NAME_COPY((list[28]).symbol, "Si");
    STR_NAME_COPY((list[28]).name, "Silicon");
    //===========================================================================
    (list[29]).proton = 14;
    (list[29]).neutron = 15;
    (list[29]).mass_of_atom = 28.976495;
    (list[29]).mass_of_nucleus = 52806.89371672695;
    (list[29]).RAM_of_element = 28.0855;
    (list[29]).abundency = 0.046832000000000006;
    STR_NAME_COPY((list[29]).symbol, "Si");
    STR_NAME_COPY((list[29]).name, "Silicon");
    //===========================================================================
    (list[30]).proton = 14;
    (list[30]).neutron = 16;
    (list[30]).mass_of_atom = 29.973770;
    (list[30]).mass_of_nucleus = 54624.813957989696;
    (list[30]).RAM_of_element = 28.0855;
    (list[30]).abundency = 0.030872;
    STR_NAME_COPY((list[30]).symbol, "Si");
    STR_NAME_COPY((list[30]).name, "Silicon");
    //===========================================================================
    (list[31]).proton = 15;
    (list[31]).neutron = 16;
    (list[31]).mass_of_atom = 30.973762;
    (list[31]).mass_of_nucleus = 56446.686984888824;
    (list[31]).RAM_of_element = 30.9738;
    (list[31]).abundency = 1.0;
    STR_NAME_COPY((list[31]).symbol, "P");
    STR_NAME_COPY((list[31]).name, "Phosphorus");
    //===========================================================================
    (list[32]).proton = 16;
    (list[32]).neutron = 16;
    (list[32]).mass_of_atom = 31.972071;
    (list[32]).mass_of_nucleus = 58265.49209194031;
    (list[32]).RAM_of_element = 32.065;
    (list[32]).abundency = 0.9493;
    STR_NAME_COPY((list[32]).symbol, "S");
    STR_NAME_COPY((list[32]).name, "Sulfur");
    //===========================================================================
    (list[33]).proton = 16;
    (list[33]).neutron = 17;
    (list[33]).mass_of_atom = 32.971458;
    (list[33]).mass_of_nucleus = 60087.262271835374;
    (list[33]).RAM_of_element = 32.065;
    (list[33]).abundency = 0.0076;
    STR_NAME_COPY((list[33]).symbol, "S");
    STR_NAME_COPY((list[33]).name, "Sulfur");
    //===========================================================================
    (list[34]).proton = 16;
    (list[34]).neutron = 18;
    (list[34]).mass_of_atom = 33.967867;
    (list[34]).mass_of_nucleus = 61903.60389242787;
    (list[34]).RAM_of_element = 32.065;
    (list[34]).abundency = 0.0429;
    STR_NAME_COPY((list[34]).symbol, "S");
    STR_NAME_COPY((list[34]).name, "Sulfur");
    //===========================================================================
    (list[35]).proton = 16;
    (list[35]).neutron = 20;
    (list[35]).mass_of_atom = 35.967081;
    (list[35]).mass_of_nucleus = 65547.94632276641;
    (list[35]).RAM_of_element = 32.065;
    (list[35]).abundency = 0.0002;
    STR_NAME_COPY((list[35]).symbol, "S");
    STR_NAME_COPY((list[35]).name, "Sulfur");
    //===========================================================================
    (list[36]).proton = 17;
    (list[36]).neutron = 18;
    (list[36]).mass_of_atom = 34.968853;
    (list[36]).mass_of_nucleus = 63727.28886961134;
    (list[36]).RAM_of_element = 35.453;
    (list[36]).abundency = 0.7578;
    STR_NAME_COPY((list[36]).symbol, "Cl");
    STR_NAME_COPY((list[36]).name, "Chlorine");
    //===========================================================================
    (list[37]).proton = 17;
    (list[37]).neutron = 20;
    (list[37]).mass_of_atom = 36.965903;
    (list[37]).mass_of_nucleus = 67367.68657116183;
    (list[37]).RAM_of_element = 35.453;
    (list[37]).abundency = 0.2422;
    STR_NAME_COPY((list[37]).symbol, "Cl");
    STR_NAME_COPY((list[37]).name, "Chlorine");
    //===========================================================================
    (list[38]).proton = 18;
    (list[38]).neutron = 18;
    (list[38]).mass_of_atom = 35.967546;
    (list[38]).mass_of_nucleus = 65546.79396550506;
    (list[38]).RAM_of_element = 39.948;
    (list[38]).abundency = 0.0033650000000000004;
    STR_NAME_COPY((list[38]).symbol, "Ar");
    STR_NAME_COPY((list[38]).name, "Argon");
    //===========================================================================
    (list[39]).proton = 18;
    (list[39]).neutron = 20;
    (list[39]).mass_of_atom = 37.962732;
    (list[39]).mass_of_nucleus = 69183.79380455053;
    (list[39]).RAM_of_element = 39.948;
    (list[39]).abundency = 0.0006320000000000001;
    STR_NAME_COPY((list[39]).symbol, "Ar");
    STR_NAME_COPY((list[39]).name, "Argon");
    //===========================================================================
    (list[40]).proton = 18;
    (list[40]).neutron = 22;
    (list[40]).mass_of_atom = 39.962383;
    (list[40]).mass_of_nucleus = 72828.93283677463;
    (list[40]).RAM_of_element = 39.948;
    (list[40]).abundency = 0.9960030000000001;
    STR_NAME_COPY((list[40]).symbol, "Ar");
    STR_NAME_COPY((list[40]).name, "Argon");
    //===========================================================================
    (list[41]).proton = 19;
    (list[41]).neutron = 20;
    (list[41]).mass_of_atom = 38.963707;
    (list[41]).mass_of_nucleus = 71007.45872997027;
    (list[41]).RAM_of_element = 39.0983;
    (list[41]).abundency = 0.932581;
    STR_NAME_COPY((list[41]).symbol, "K");
    STR_NAME_COPY((list[41]).name, "Potassium");
    //===========================================================================
    (list[42]).proton = 19;
    (list[42]).neutron = 21;
    (list[42]).mass_of_atom = 39.963999;
    (list[42]).mass_of_nucleus = 72830.8786231524;
    (list[42]).RAM_of_element = 39.0983;
    (list[42]).abundency = 0.000117;
    STR_NAME_COPY((list[42]).symbol, "K");
    STR_NAME_COPY((list[42]).name, "Potassium");
    //===========================================================================
    (list[43]).proton = 19;
    (list[43]).neutron = 22;
    (list[43]).mass_of_atom = 40.961826;
    (list[43]).mass_of_nucleus = 74649.80509837586;
    (list[43]).RAM_of_element = 39.0983;
    (list[43]).abundency = 0.067302;
    STR_NAME_COPY((list[43]).symbol, "K");
    STR_NAME_COPY((list[43]).name, "Potassium");
    //===========================================================================
    (list[44]).proton = 20;
    (list[44]).neutron = 20;
    (list[44]).mass_of_atom = 39.962591;
    (list[44]).mass_of_nucleus = 72827.31199739751;
    (list[44]).RAM_of_element = 40.078;
    (list[44]).abundency = 0.96941;
    STR_NAME_COPY((list[44]).symbol, "Ca");
    STR_NAME_COPY((list[44]).name, "Calcium");
    //===========================================================================
    (list[45]).proton = 20;
    (list[45]).neutron = 22;
    (list[45]).mass_of_atom = 41.958618;
    (list[45]).mass_of_nucleus = 76465.84488492298;
    (list[45]).RAM_of_element = 40.078;
    (list[45]).abundency = 0.00647;
    STR_NAME_COPY((list[45]).symbol, "Ca");
    STR_NAME_COPY((list[45]).name, "Calcium");
    //===========================================================================
    (list[46]).proton = 20;
    (list[46]).neutron = 23;
    (list[46]).mass_of_atom = 42.958767;
    (list[46]).mass_of_nucleus = 78289.00410517688;
    (list[46]).RAM_of_element = 40.078;
    (list[46]).abundency = 0.00135;
    STR_NAME_COPY((list[46]).symbol, "Ca");
    STR_NAME_COPY((list[46]).name, "Calcium");
    //===========================================================================
    (list[47]).proton = 20;
    (list[47]).neutron = 24;
    (list[47]).mass_of_atom = 43.955481;
    (list[47]).mass_of_nucleus = 80105.9017064904;
    (list[47]).RAM_of_element = 40.078;
    (list[47]).abundency = 0.02086;
    STR_NAME_COPY((list[47]).symbol, "Ca");
    STR_NAME_COPY((list[47]).name, "Calcium");
    //===========================================================================
    (list[48]).proton = 20;
    (list[48]).neutron = 26;
    (list[48]).mass_of_atom = 45.953693;
    (list[48]).mass_of_nucleus = 83748.41760344373;
    (list[48]).RAM_of_element = 40.078;
    (list[48]).abundency = 4e-05;
    STR_NAME_COPY((list[48]).symbol, "Ca");
    STR_NAME_COPY((list[48]).name, "Calcium");
    //===========================================================================
    (list[49]).proton = 20;
    (list[49]).neutron = 28;
    (list[49]).mass_of_atom = 47.952534;
    (list[49]).mass_of_nucleus = 87392.08009670374;
    (list[49]).RAM_of_element = 40.078;
    (list[49]).abundency = 0.00187;
    STR_NAME_COPY((list[49]).symbol, "Ca");
    STR_NAME_COPY((list[49]).name, "Calcium");
    //===========================================================================
    (list[50]).proton = 21;
    (list[50]).neutron = 24;
    (list[50]).mass_of_atom = 44.955910;
    (list[50]).mass_of_nucleus = 81928.5713352751;
    (list[50]).RAM_of_element = 44.9559;
    (list[50]).abundency = 1.0;
    STR_NAME_COPY((list[50]).symbol, "Sc");
    STR_NAME_COPY((list[50]).name, "Scandium");
    //===========================================================================
    (list[51]).proton = 22;
    (list[51]).neutron = 24;
    (list[51]).mass_of_atom = 45.952629;
    (list[51]).mass_of_nucleus = 83744.47805102669;
    (list[51]).RAM_of_element = 47.867;
    (list[51]).abundency = 0.0825;
    STR_NAME_COPY((list[51]).symbol, "Ti");
    STR_NAME_COPY((list[51]).name, "Titanium");
    //===========================================================================
    (list[52]).proton = 22;
    (list[52]).neutron = 25;
    (list[52]).mass_of_atom = 46.951764;
    (list[52]).mass_of_nucleus = 85565.78886324403;
    (list[52]).RAM_of_element = 47.867;
    (list[52]).abundency = 0.07440000000000001;
    STR_NAME_COPY((list[52]).symbol, "Ti");
    STR_NAME_COPY((list[52]).name, "Titanium");
    //===========================================================================
    (list[53]).proton = 22;
    (list[53]).neutron = 26;
    (list[53]).mass_of_atom = 47.947947;
    (list[53]).mass_of_nucleus = 87381.71851123667;
    (list[53]).RAM_of_element = 47.867;
    (list[53]).abundency = 0.7372;
    STR_NAME_COPY((list[53]).symbol, "Ti");
    STR_NAME_COPY((list[53]).name, "Titanium");
    //===========================================================================
    (list[54]).proton = 22;
    (list[54]).neutron = 27;
    (list[54]).mass_of_atom = 48.947871;
    (list[54]).mass_of_nucleus = 89204.46758177831;
    (list[54]).RAM_of_element = 47.867;
    (list[54]).abundency = 0.0541;
    STR_NAME_COPY((list[54]).symbol, "Ti");
    STR_NAME_COPY((list[54]).name, "Titanium");
    //===========================================================================
    (list[55]).proton = 22;
    (list[55]).neutron = 28;
    (list[55]).mass_of_atom = 49.944792;
    (list[55]).mass_of_nucleus = 91021.74252082712;
    (list[55]).RAM_of_element = 47.867;
    (list[55]).abundency = 0.0518;
    STR_NAME_COPY((list[55]).symbol, "Ti");
    STR_NAME_COPY((list[55]).name, "Titanium");
    //===========================================================================
    (list[56]).proton = 23;
    (list[56]).neutron = 27;
    (list[56]).mass_of_atom = 49.947163;
    (list[56]).mass_of_nucleus = 91025.06458735044;
    (list[56]).RAM_of_element = 50.9415;
    (list[56]).abundency = 0.0025;
    STR_NAME_COPY((list[56]).symbol, "V");
    STR_NAME_COPY((list[56]).name, "Vanadium");
    //===========================================================================
    (list[57]).proton = 23;
    (list[57]).neutron = 28;
    (list[57]).mass_of_atom = 50.943964;
    (list[57]).mass_of_nucleus = 92842.12077988604;
    (list[57]).RAM_of_element = 50.9415;
    (list[57]).abundency = 0.9975;
    STR_NAME_COPY((list[57]).symbol, "V");
    STR_NAME_COPY((list[57]).name, "Vanadium");
    //===========================================================================
    (list[58]).proton = 24;
    (list[58]).neutron = 26;
    (list[58]).mass_of_atom = 49.946050;
    (list[58]).mass_of_nucleus = 91022.0357134405;
    (list[58]).RAM_of_element = 51.9961;
    (list[58]).abundency = 0.043449999999999996;
    STR_NAME_COPY((list[58]).symbol, "Cr");
    STR_NAME_COPY((list[58]).name, "Chromium");
    //===========================================================================
    (list[59]).proton = 24;
    (list[59]).neutron = 28;
    (list[59]).mass_of_atom = 51.940512;
    (list[59]).mass_of_nucleus = 94657.71578185631;
    (list[59]).RAM_of_element = 51.9961;
    (list[59]).abundency = 0.83789;
    STR_NAME_COPY((list[59]).symbol, "Cr");
    STR_NAME_COPY((list[59]).name, "Chromium");
    //===========================================================================
    (list[60]).proton = 24;
    (list[60]).neutron = 29;
    (list[60]).mass_of_atom = 52.940654;
    (list[60]).mass_of_nucleus = 96480.86224189695;
    (list[60]).RAM_of_element = 51.9961;
    (list[60]).abundency = 0.09501;
    STR_NAME_COPY((list[60]).symbol, "Cr");
    STR_NAME_COPY((list[60]).name, "Chromium");
    //===========================================================================
    (list[61]).proton = 24;
    (list[61]).neutron = 30;
    (list[61]).mass_of_atom = 53.938885;
    (list[61]).mass_of_nucleus = 98300.52516371485;
    (list[61]).RAM_of_element = 51.9961;
    (list[61]).abundency = 0.02365;
    STR_NAME_COPY((list[61]).symbol, "Cr");
    STR_NAME_COPY((list[61]).name, "Chromium");
    //===========================================================================
    (list[62]).proton = 25;
    (list[62]).neutron = 30;
    (list[62]).mass_of_atom = 54.938050;
    (list[62]).mass_of_nucleus = 100120.8906625605;
    (list[62]).RAM_of_element = 54.938;
    (list[62]).abundency = 1.0;
    STR_NAME_COPY((list[62]).symbol, "Mn");
    STR_NAME_COPY((list[62]).name, "Manganese");
    //===========================================================================
    (list[63]).proton = 26;
    (list[63]).neutron = 28;
    (list[63]).mass_of_atom = 53.939615;
    (list[63]).mass_of_nucleus = 98299.85587167015;
    (list[63]).RAM_of_element = 55.845;
    (list[63]).abundency = 0.058449999999999995;
    STR_NAME_COPY((list[63]).symbol, "Fe");
    STR_NAME_COPY((list[63]).name, "Iron");
    //===========================================================================
    (list[64]).proton = 26;
    (list[64]).neutron = 30;
    (list[64]).mass_of_atom = 55.934942;
    (list[64]).mass_of_nucleus = 101937.11273786862;
    (list[64]).RAM_of_element = 55.845;
    (list[64]).abundency = 0.91754;
    STR_NAME_COPY((list[64]).symbol, "Fe");
    STR_NAME_COPY((list[64]).name, "Iron");
    //===========================================================================
    (list[65]).proton = 26;
    (list[65]).neutron = 31;
    (list[65]).mass_of_atom = 56.935399;
    (list[65]).mass_of_nucleus = 103760.83340750639;
    (list[65]).RAM_of_element = 55.845;
    (list[65]).abundency = 0.02119;
    STR_NAME_COPY((list[65]).symbol, "Fe");
    STR_NAME_COPY((list[65]).name, "Iron");
    //===========================================================================
    (list[66]).proton = 26;
    (list[66]).neutron = 32;
    (list[66]).mass_of_atom = 57.933280;
    (list[66]).mass_of_nucleus = 105579.8583186608;
    (list[66]).RAM_of_element = 55.845;
    (list[66]).abundency = 0.0028199999999999996;
    STR_NAME_COPY((list[66]).symbol, "Fe");
    STR_NAME_COPY((list[66]).name, "Iron");
    //===========================================================================
    (list[67]).proton = 27;
    (list[67]).neutron = 32;
    (list[67]).mass_of_atom = 58.933200;
    (list[67]).mass_of_nucleus = 107401.600097652;
    (list[67]).RAM_of_element = 58.9332;
    (list[67]).abundency = 1.0;
    STR_NAME_COPY((list[67]).symbol, "Co");
    STR_NAME_COPY((list[67]).name, "Cobalt");
    //===========================================================================
    (list[68]).proton = 28;
    (list[68]).neutron = 30;
    (list[68]).mass_of_atom = 57.935348;
    (list[68]).mass_of_nucleus = 105581.62805023827;
    (list[68]).RAM_of_element = 58.6934;
    (list[68]).abundency = 0.680769;
    STR_NAME_COPY((list[68]).symbol, "Ni");
    STR_NAME_COPY((list[68]).name, "Nickel");
    //===========================================================================
    (list[69]).proton = 28;
    (list[69]).neutron = 32;
    (list[69]).mass_of_atom = 59.930791;
    (list[69]).mass_of_nucleus = 109219.09637139952;
    (list[69]).RAM_of_element = 58.6934;
    (list[69]).abundency = 0.262231;
    STR_NAME_COPY((list[69]).symbol, "Ni");
    STR_NAME_COPY((list[69]).name, "Nickel");
    //===========================================================================
    (list[70]).proton = 28;
    (list[70]).neutron = 33;
    (list[70]).mass_of_atom = 60.931060;
    (list[70]).mass_of_nucleus = 111042.4743381666;
    (list[70]).RAM_of_element = 58.6934;
    (list[70]).abundency = 0.011399;
    STR_NAME_COPY((list[70]).symbol, "Ni");
    STR_NAME_COPY((list[70]).name, "Nickel");
    //===========================================================================
    (list[71]).proton = 28;
    (list[71]).neutron = 34;
    (list[71]).mass_of_atom = 61.928349;
    (list[71]).mass_of_nucleus = 112860.42009985588;
    (list[71]).RAM_of_element = 58.6934;
    (list[71]).abundency = 0.036345;
    STR_NAME_COPY((list[71]).symbol, "Ni");
    STR_NAME_COPY((list[71]).name, "Nickel");
    //===========================================================================
    (list[72]).proton = 28;
    (list[72]).neutron = 36;
    (list[72]).mass_of_atom = 63.927970;
    (list[72]).mass_of_nucleus = 116505.5044454517;
    (list[72]).RAM_of_element = 58.6934;
    (list[72]).abundency = 0.009256;
    STR_NAME_COPY((list[72]).symbol, "Ni");
    STR_NAME_COPY((list[72]).name, "Nickel");
    //===========================================================================
    (list[73]).proton = 29;
    (list[73]).neutron = 34;
    (list[73]).mass_of_atom = 62.929601;
    (list[73]).mass_of_nucleus = 114684.5899651436;
    (list[73]).RAM_of_element = 63.546;
    (list[73]).abundency = 0.6917;
    STR_NAME_COPY((list[73]).symbol, "Cu");
    STR_NAME_COPY((list[73]).name, "Copper");
    //===========================================================================
    (list[74]).proton = 29;
    (list[74]).neutron = 36;
    (list[74]).mass_of_atom = 64.927794;
    (list[74]).mass_of_nucleus = 118327.07122723234;
    (list[74]).RAM_of_element = 63.546;
    (list[74]).abundency = 0.30829999999999996;
    STR_NAME_COPY((list[74]).symbol, "Cu");
    STR_NAME_COPY((list[74]).name, "Copper");
    //===========================================================================
    (list[75]).proton = 30;
    (list[75]).neutron = 34;
    (list[75]).mass_of_atom = 63.929147;
    (list[75]).mass_of_nucleus = 116505.64998416867;
    (list[75]).RAM_of_element = 65.39;
    (list[75]).abundency = 0.4863;
    STR_NAME_COPY((list[75]).symbol, "Zn");
    STR_NAME_COPY((list[75]).name, "Zinc");
    //===========================================================================
    (list[76]).proton = 30;
    (list[76]).neutron = 36;
    (list[76]).mass_of_atom = 65.926037;
    (list[76]).mass_of_nucleus = 120145.75602370156;
    (list[76]).RAM_of_element = 65.39;
    (list[76]).abundency = 0.27899999999999997;
    STR_NAME_COPY((list[76]).symbol, "Zn");
    STR_NAME_COPY((list[76]).name, "Zinc");
    //===========================================================================
    (list[77]).proton = 30;
    (list[77]).neutron = 37;
    (list[77]).mass_of_atom = 66.927131;
    (list[77]).mass_of_nucleus = 121970.63787274691;
    (list[77]).RAM_of_element = 65.39;
    (list[77]).abundency = 0.040999999999999995;
    STR_NAME_COPY((list[77]).symbol, "Zn");
    STR_NAME_COPY((list[77]).name, "Zinc");
    //===========================================================================
    (list[78]).proton = 30;
    (list[78]).neutron = 38;
    (list[78]).mass_of_atom = 67.924848;
    (list[78]).mass_of_nucleus = 123789.36383033327;
    (list[78]).RAM_of_element = 65.39;
    (list[78]).abundency = 0.1875;
    STR_NAME_COPY((list[78]).symbol, "Zn");
    STR_NAME_COPY((list[78]).name, "Zinc");
    //===========================================================================
    (list[79]).proton = 30;
    (list[79]).neutron = 40;
    (list[79]).mass_of_atom = 69.925325;
    (list[79]).mass_of_nucleus = 127436.00856772326;
    (list[79]).RAM_of_element = 65.39;
    (list[79]).abundency = 0.0062;
    STR_NAME_COPY((list[79]).symbol, "Zn");
    STR_NAME_COPY((list[79]).name, "Zinc");
    //===========================================================================
    (list[80]).proton = 31;
    (list[80]).neutron = 38;
    (list[80]).mass_of_atom = 68.925581;
    (list[80]).mass_of_nucleus = 125612.5876169514;
    (list[80]).RAM_of_element = 69.723;
    (list[80]).abundency = 0.60108;
    STR_NAME_COPY((list[80]).symbol, "Ga");
    STR_NAME_COPY((list[80]).name, "Gallium");
    //===========================================================================
    (list[81]).proton = 31;
    (list[81]).neutron = 40;
    (list[81]).mass_of_atom = 70.924705;
    (list[81]).mass_of_nucleus = 129256.76598740506;
    (list[81]).RAM_of_element = 69.723;
    (list[81]).abundency = 0.39892000000000005;
    STR_NAME_COPY((list[81]).symbol, "Ga");
    STR_NAME_COPY((list[81]).name, "Gallium");
    //===========================================================================
    (list[82]).proton = 32;
    (list[82]).neutron = 38;
    (list[82]).mass_of_atom = 69.924250;
    (list[82]).mass_of_nucleus = 127432.0489635425;
    (list[82]).RAM_of_element = 72.64;
    (list[82]).abundency = 0.2084;
    STR_NAME_COPY((list[82]).symbol, "Ge");
    STR_NAME_COPY((list[82]).name, "Germanium");
    //===========================================================================
    (list[83]).proton = 32;
    (list[83]).neutron = 40;
    (list[83]).mass_of_atom = 71.922076;
    (list[83]).mass_of_nucleus = 131073.86122587835;
    (list[83]).RAM_of_element = 72.64;
    (list[83]).abundency = 0.2754;
    STR_NAME_COPY((list[83]).symbol, "Ge");
    STR_NAME_COPY((list[83]).name, "Germanium");
    //===========================================================================
    (list[84]).proton = 32;
    (list[84]).neutron = 41;
    (list[84]).mass_of_atom = 72.923459;
    (list[84]).mass_of_nucleus = 132899.26988944298;
    (list[84]).RAM_of_element = 72.64;
    (list[84]).abundency = 0.07730000000000001;
    STR_NAME_COPY((list[84]).symbol, "Ge");
    STR_NAME_COPY((list[84]).name, "Germanium");
    //===========================================================================
    (list[85]).proton = 32;
    (list[85]).neutron = 42;
    (list[85]).mass_of_atom = 73.921178;
    (list[85]).mass_of_nucleus = 134717.99949280458;
    (list[85]).RAM_of_element = 72.64;
    (list[85]).abundency = 0.3628;
    STR_NAME_COPY((list[85]).symbol, "Ge");
    STR_NAME_COPY((list[85]).name, "Germanium");
    //===========================================================================
    (list[86]).proton = 32;
    (list[86]).neutron = 44;
    (list[86]).mass_of_atom = 75.921403;
    (list[86]).mass_of_nucleus = 138364.18486251682;
    (list[86]).RAM_of_element = 72.64;
    (list[86]).abundency = 0.0761;
    STR_NAME_COPY((list[86]).symbol, "Ge");
    STR_NAME_COPY((list[86]).name, "Germanium");
    //===========================================================================
    (list[87]).proton = 33;
    (list[87]).neutron = 42;
    (list[87]).mass_of_atom = 74.921596;
    (list[87]).mass_of_nucleus = 136540.64906982554;
    (list[87]).RAM_of_element = 74.9216;
    (list[87]).abundency = 1.0;
    STR_NAME_COPY((list[87]).symbol, "As");
    STR_NAME_COPY((list[87]).name, "Arsenic");
    //===========================================================================
    (list[88]).proton = 34;
    (list[88]).neutron = 40;
    (list[88]).mass_of_atom = 73.922477;
    (list[88]).mass_of_nucleus = 134718.36742380998;
    (list[88]).RAM_of_element = 78.96;
    (list[88]).abundency = 0.0089;
    STR_NAME_COPY((list[88]).symbol, "Se");
    STR_NAME_COPY((list[88]).name, "Selenium");
    //===========================================================================
    (list[89]).proton = 34;
    (list[89]).neutron = 42;
    (list[89]).mass_of_atom = 75.919214;
    (list[89]).mass_of_nucleus = 138358.19456153852;
    (list[89]).RAM_of_element = 78.96;
    (list[89]).abundency = 0.09369999999999999;
    STR_NAME_COPY((list[89]).symbol, "Se");
    STR_NAME_COPY((list[89]).name, "Selenium");
    //===========================================================================
    (list[90]).proton = 34;
    (list[90]).neutron = 43;
    (list[90]).mass_of_atom = 76.919915;
    (list[90]).mass_of_nucleus = 140182.36001575316;
    (list[90]).RAM_of_element = 78.96;
    (list[90]).abundency = 0.07629999999999999;
    STR_NAME_COPY((list[90]).symbol, "Se");
    STR_NAME_COPY((list[90]).name, "Selenium");
    //===========================================================================
    (list[91]).proton = 34;
    (list[91]).neutron = 44;
    (list[91]).mass_of_atom = 77.917310;
    (list[91]).mass_of_nucleus = 142000.4990035291;
    (list[91]).RAM_of_element = 78.96;
    (list[91]).abundency = 0.2377;
    STR_NAME_COPY((list[91]).symbol, "Se");
    STR_NAME_COPY((list[91]).name, "Selenium");
    //===========================================================================
    (list[92]).proton = 34;
    (list[92]).neutron = 46;
    (list[92]).mass_of_atom = 79.916522;
    (list[92]).mass_of_nucleus = 145644.8377880924;
    (list[92]).RAM_of_element = 78.96;
    (list[92]).abundency = 0.4961;
    STR_NAME_COPY((list[92]).symbol, "Se");
    STR_NAME_COPY((list[92]).name, "Selenium");
    //===========================================================================
    (list[93]).proton = 34;
    (list[93]).neutron = 48;
    (list[93]).mass_of_atom = 81.916700;
    (list[93]).mass_of_nucleus = 149290.937482087;
    (list[93]).RAM_of_element = 78.96;
    (list[93]).abundency = 0.0873;
    STR_NAME_COPY((list[93]).symbol, "Se");
    STR_NAME_COPY((list[93]).name, "Selenium");
    //===========================================================================
    (list[94]).proton = 35;
    (list[94]).neutron = 44;
    (list[94]).mass_of_atom = 78.918338;
    (list[94]).mass_of_nucleus = 143824.2605419922;
    (list[94]).RAM_of_element = 79.904;
    (list[94]).abundency = 0.5069;
    STR_NAME_COPY((list[94]).symbol, "Br");
    STR_NAME_COPY((list[94]).name, "Bromine");
    //===========================================================================
    (list[95]).proton = 35;
    (list[95]).neutron = 46;
    (list[95]).mass_of_atom = 80.916291;
    (list[95]).mass_of_nucleus = 147466.30431105453;
    (list[95]).RAM_of_element = 79.904;
    (list[95]).abundency = 0.49310000000000004;
    STR_NAME_COPY((list[95]).symbol, "Br");
    STR_NAME_COPY((list[95]).name, "Bromine");
    //===========================================================================
    (list[96]).proton = 36;
    (list[96]).neutron = 42;
    (list[96]).mass_of_atom = 77.920386;
    (list[96]).mass_of_nucleus = 142004.10620581746;
    (list[96]).RAM_of_element = 83.8;
    (list[96]).abundency = 0.0034999999999999996;
    STR_NAME_COPY((list[96]).symbol, "Kr");
    STR_NAME_COPY((list[96]).name, "Krypton");
    //===========================================================================
    (list[97]).proton = 36;
    (list[97]).neutron = 44;
    (list[97]).mass_of_atom = 79.916378;
    (list[97]).mass_of_nucleus = 145642.57529227657;
    (list[97]).RAM_of_element = 83.8;
    (list[97]).abundency = 0.022799999999999997;
    STR_NAME_COPY((list[97]).symbol, "Kr");
    STR_NAME_COPY((list[97]).name, "Krypton");
    //===========================================================================
    (list[98]).proton = 36;
    (list[98]).neutron = 46;
    (list[98]).mass_of_atom = 81.913485;
    (list[98]).mass_of_nucleus = 149283.07689842084;
    (list[98]).RAM_of_element = 83.8;
    (list[98]).abundency = 0.1158;
    STR_NAME_COPY((list[98]).symbol, "Kr");
    STR_NAME_COPY((list[98]).name, "Krypton");
    //===========================================================================
    (list[99]).proton = 36;
    (list[99]).neutron = 47;
    (list[99]).mass_of_atom = 82.914136;
    (list[99]).mass_of_nucleus = 151107.15120825495;
    (list[99]).RAM_of_element = 83.8;
    (list[99]).abundency = 0.1149;
    STR_NAME_COPY((list[99]).symbol, "Kr");
    STR_NAME_COPY((list[99]).name, "Krypton");
    //===========================================================================
    (list[100]).proton = 36;
    (list[100]).neutron = 48;
    (list[100]).mass_of_atom = 83.911507;
    (list[100]).mass_of_nucleus = 152925.24644672827;
    (list[100]).RAM_of_element = 83.8;
    (list[100]).abundency = 0.57;
    STR_NAME_COPY((list[100]).symbol, "Kr");
    STR_NAME_COPY((list[100]).name, "Krypton");
    //===========================================================================
    (list[101]).proton = 36;
    (list[101]).neutron = 50;
    (list[101]).mass_of_atom = 85.910610;
    (list[101]).mass_of_nucleus = 156569.38653654212;
    (list[101]).RAM_of_element = 83.8;
    (list[101]).abundency = 0.17300000000000001;
    STR_NAME_COPY((list[101]).symbol, "Kr");
    STR_NAME_COPY((list[101]).name, "Krypton");
    //===========================================================================
    (list[102]).proton = 37;
    (list[102]).neutron = 48;
    (list[102]).mass_of_atom = 84.911789;
    (list[102]).mass_of_nucleus = 154747.6481110343;
    (list[102]).RAM_of_element = 85.4678;
    (list[102]).abundency = 0.7217;
    STR_NAME_COPY((list[102]).symbol, "Rb");
    STR_NAME_COPY((list[102]).name, "Rubidium");
    //===========================================================================
    (list[103]).proton = 37;
    (list[103]).neutron = 50;
    (list[103]).mass_of_atom = 86.909183;
    (list[103]).mass_of_nucleus = 158388.67288592263;
    (list[103]).RAM_of_element = 85.4678;
    (list[103]).abundency = 0.2783;
    STR_NAME_COPY((list[103]).symbol, "Rb");
    STR_NAME_COPY((list[103]).name, "Rubidium");
    //===========================================================================
    (list[104]).proton = 38;
    (list[104]).neutron = 46;
    (list[104]).mass_of_atom = 83.913425;
    (list[104]).mass_of_nucleus = 152926.74274516426;
    (list[104]).RAM_of_element = 87.62;
    (list[104]).abundency = 0.005600000000000001;
    STR_NAME_COPY((list[104]).symbol, "Sr");
    STR_NAME_COPY((list[104]).name, "Strontium");
    //===========================================================================
    (list[105]).proton = 38;
    (list[105]).neutron = 48;
    (list[105]).mass_of_atom = 85.909262;
    (list[105]).mass_of_nucleus = 156564.9292840438;
    (list[105]).RAM_of_element = 87.62;
    (list[105]).abundency = 0.0986;
    STR_NAME_COPY((list[105]).symbol, "Sr");
    STR_NAME_COPY((list[105]).name, "Strontium");
    //===========================================================================
    (list[106]).proton = 38;
    (list[106]).neutron = 49;
    (list[106]).mass_of_atom = 86.908879;
    (list[106]).mass_of_nucleus = 158387.11872808918;
    (list[106]).RAM_of_element = 87.62;
    (list[106]).abundency = 0.07;
    STR_NAME_COPY((list[106]).symbol, "Sr");
    STR_NAME_COPY((list[106]).name, "Strontium");
    //===========================================================================
    (list[107]).proton = 38;
    (list[107]).neutron = 50;
    (list[107]).mass_of_atom = 87.905614;
    (list[107]).mass_of_nucleus = 160204.05461004254;
    (list[107]).RAM_of_element = 87.62;
    (list[107]).abundency = 0.8258;
    STR_NAME_COPY((list[107]).symbol, "Sr");
    STR_NAME_COPY((list[107]).name, "Strontium");
    //===========================================================================
    (list[108]).proton = 39;
    (list[108]).neutron = 50;
    (list[108]).mass_of_atom = 88.905848;
    (list[108]).mass_of_nucleus = 162026.3687757433;
    (list[108]).RAM_of_element = 88.9059;
    (list[108]).abundency = 1.0;
    STR_NAME_COPY((list[108]).symbol, "Y");
    STR_NAME_COPY((list[108]).name, "Yttrium");
    //===========================================================================
    (list[109]).proton = 40;
    (list[109]).neutron = 50;
    (list[109]).mass_of_atom = 89.904704;
    (list[109]).mass_of_nucleus = 163846.17100231742;
    (list[109]).RAM_of_element = 91.224;
    (list[109]).abundency = 0.5145000000000001;
    STR_NAME_COPY((list[109]).symbol, "Zr");
    STR_NAME_COPY((list[109]).name, "Zirconium");
    //===========================================================================
    (list[110]).proton = 40;
    (list[110]).neutron = 51;
    (list[110]).mass_of_atom = 90.905645;
    (list[110]).mass_of_nucleus = 165670.77394955847;
    (list[110]).RAM_of_element = 91.224;
    (list[110]).abundency = 0.11220000000000001;
    STR_NAME_COPY((list[110]).symbol, "Zr");
    STR_NAME_COPY((list[110]).name, "Zirconium");
    //===========================================================================
    (list[111]).proton = 40;
    (list[111]).neutron = 52;
    (list[111]).mass_of_atom = 91.905040;
    (list[111]).mass_of_nucleus = 167492.5587125544;
    (list[111]).RAM_of_element = 91.224;
    (list[111]).abundency = 0.17149999999999999;
    STR_NAME_COPY((list[111]).symbol, "Zr");
    STR_NAME_COPY((list[111]).name, "Zirconium");
    //===========================================================================
    (list[112]).proton = 40;
    (list[112]).neutron = 54;
    (list[112]).mass_of_atom = 93.906316;
    (list[112]).mass_of_nucleus = 171140.65993714478;
    (list[112]).RAM_of_element = 91.224;
    (list[112]).abundency = 0.17379999999999998;
    STR_NAME_COPY((list[112]).symbol, "Zr");
    STR_NAME_COPY((list[112]).name, "Zirconium");
    //===========================================================================
    (list[113]).proton = 40;
    (list[113]).neutron = 56;
    (list[113]).mass_of_atom = 95.908276;
    (list[113]).mass_of_nucleus = 174790.00801686037;
    (list[113]).RAM_of_element = 91.224;
    (list[113]).abundency = 0.027999999999999997;
    STR_NAME_COPY((list[113]).symbol, "Zr");
    STR_NAME_COPY((list[113]).name, "Zirconium");
    //===========================================================================
    (list[114]).proton = 41;
    (list[114]).neutron = 52;
    (list[114]).mass_of_atom = 92.906378;
    (list[114]).mass_of_nucleus = 169316.8853461766;
    (list[114]).RAM_of_element = 92.9064;
    (list[114]).abundency = 1.0;
    STR_NAME_COPY((list[114]).symbol, "Nb");
    STR_NAME_COPY((list[114]).name, "Niobium");
    //===========================================================================
    (list[115]).proton = 42;
    (list[115]).neutron = 50;
    (list[115]).mass_of_atom = 91.906810;
    (list[115]).mass_of_nucleus = 167493.7852236241;
    (list[115]).RAM_of_element = 95.94;
    (list[115]).abundency = 0.1484;
    STR_NAME_COPY((list[115]).symbol, "Mo");
    STR_NAME_COPY((list[115]).name, "Molybdenum");
    //===========================================================================
    (list[116]).proton = 42;
    (list[116]).neutron = 52;
    (list[116]).mass_of_atom = 93.905088;
    (list[116]).mass_of_nucleus = 171136.4214311597;
    (list[116]).RAM_of_element = 95.94;
    (list[116]).abundency = 0.0925;
    STR_NAME_COPY((list[116]).symbol, "Mo");
    STR_NAME_COPY((list[116]).name, "Molybdenum");
    //===========================================================================
    (list[117]).proton = 42;
    (list[117]).neutron = 53;
    (list[117]).mass_of_atom = 94.905841;
    (list[117]).mass_of_nucleus = 172960.68167553;
    (list[117]).RAM_of_element = 95.94;
    (list[117]).abundency = 0.1592;
    STR_NAME_COPY((list[117]).symbol, "Mo");
    STR_NAME_COPY((list[117]).name, "Molybdenum");
    //===========================================================================
    (list[118]).proton = 42;
    (list[118]).neutron = 54;
    (list[118]).mass_of_atom = 95.904679;
    (list[118]).mass_of_nucleus = 174781.4510901272;
    (list[118]).RAM_of_element = 95.94;
    (list[118]).abundency = 0.1668;
    STR_NAME_COPY((list[118]).symbol, "Mo");
    STR_NAME_COPY((list[118]).name, "Molybdenum");
    //===========================================================================
    (list[119]).proton = 42;
    (list[119]).neutron = 55;
    (list[119]).mass_of_atom = 96.906021;
    (list[119]).mass_of_nucleus = 176606.7850152998;
    (list[119]).RAM_of_element = 95.94;
    (list[119]).abundency = 0.0955;
    STR_NAME_COPY((list[119]).symbol, "Mo");
    STR_NAME_COPY((list[119]).name, "Molybdenum");
    //===========================================================================
    (list[120]).proton = 42;
    (list[120]).neutron = 56;
    (list[120]).mass_of_atom = 97.905408;
    (list[120]).mass_of_nucleus = 178428.55519519487;
    (list[120]).RAM_of_element = 95.94;
    (list[120]).abundency = 0.2413;
    STR_NAME_COPY((list[120]).symbol, "Mo");
    STR_NAME_COPY((list[120]).name, "Molybdenum");
    //===========================================================================
    (list[121]).proton = 42;
    (list[121]).neutron = 58;
    (list[121]).mass_of_atom = 99.907477;
    (list[121]).mass_of_nucleus = 182078.10196965997;
    (list[121]).RAM_of_element = 95.94;
    (list[121]).abundency = 0.09630000000000001;
    STR_NAME_COPY((list[121]).symbol, "Mo");
    STR_NAME_COPY((list[121]).name, "Molybdenum");
    //===========================================================================
    (list[122]).proton = 43;
    (list[122]).neutron = 55;
    (list[122]).mass_of_atom = 97.907216;
    (list[122]).mass_of_nucleus = 178430.85097599376;
    (list[122]).RAM_of_element = 98;
    (list[122]).abundency = 0.0;
    STR_NAME_COPY((list[122]).symbol, "Tc");
    STR_NAME_COPY((list[122]).name, "Technetium");
    //===========================================================================
    (list[123]).proton = 44;
    (list[123]).neutron = 52;
    (list[123]).mass_of_atom = 95.907598;
    (list[123]).mass_of_nucleus = 174784.77209906076;
    (list[123]).RAM_of_element = 101.07;
    (list[123]).abundency = 0.0554;
    STR_NAME_COPY((list[123]).symbol, "Ru");
    STR_NAME_COPY((list[123]).name, "Ruthenium");
    //===========================================================================
    (list[124]).proton = 44;
    (list[124]).neutron = 54;
    (list[124]).mass_of_atom = 97.905287;
    (list[124]).mass_of_nucleus = 178426.33462579406;
    (list[124]).RAM_of_element = 101.07;
    (list[124]).abundency = 0.0187;
    STR_NAME_COPY((list[124]).symbol, "Ru");
    STR_NAME_COPY((list[124]).name, "Ruthenium");
    //===========================================================================
    (list[125]).proton = 44;
    (list[125]).neutron = 55;
    (list[125]).mass_of_atom = 98.905939;
    (list[125]).mass_of_nucleus = 180250.4107585158;
    (list[125]).RAM_of_element = 101.07;
    (list[125]).abundency = 0.1276;
    STR_NAME_COPY((list[125]).symbol, "Ru");
    STR_NAME_COPY((list[125]).name, "Ruthenium");
    //===========================================================================
    (list[126]).proton = 44;
    (list[126]).neutron = 56;
    (list[126]).mass_of_atom = 99.904220;
    (list[126]).mass_of_nucleus = 182070.1648247142;
    (list[126]).RAM_of_element = 101.07;
    (list[126]).abundency = 0.126;
    STR_NAME_COPY((list[126]).symbol, "Ru");
    STR_NAME_COPY((list[126]).name, "Ruthenium");
    //===========================================================================
    (list[127]).proton = 44;
    (list[127]).neutron = 57;
    (list[127]).mass_of_atom = 100.905582;
    (list[127]).mass_of_nucleus = 183895.535207639;
    (list[127]).RAM_of_element = 101.07;
    (list[127]).abundency = 0.17059999999999997;
    STR_NAME_COPY((list[127]).symbol, "Ru");
    STR_NAME_COPY((list[127]).name, "Ruthenium");
    //===========================================================================
    (list[128]).proton = 44;
    (list[128]).neutron = 58;
    (list[128]).mass_of_atom = 101.904350;
    (list[128]).mass_of_nucleus = 185716.1770201035;
    (list[128]).RAM_of_element = 101.07;
    (list[128]).abundency = 0.3155;
    STR_NAME_COPY((list[128]).symbol, "Ru");
    STR_NAME_COPY((list[128]).name, "Ruthenium");
    //===========================================================================
    (list[129]).proton = 44;
    (list[129]).neutron = 60;
    (list[129]).mass_of_atom = 103.905430;
    (list[129]).mass_of_nucleus = 189363.9209587223;
    (list[129]).RAM_of_element = 101.07;
    (list[129]).abundency = 0.1862;
    STR_NAME_COPY((list[129]).symbol, "Ru");
    STR_NAME_COPY((list[129]).name, "Ruthenium");
    //===========================================================================
    (list[130]).proton = 45;
    (list[130]).neutron = 58;
    (list[130]).mass_of_atom = 102.905504;
    (list[130]).mass_of_nucleus = 187540.16824240543;
    (list[130]).RAM_of_element = 102.9055;
    (list[130]).abundency = 1.0;
    STR_NAME_COPY((list[130]).symbol, "Rh");
    STR_NAME_COPY((list[130]).name, "Rhodium");
    //===========================================================================
    (list[131]).proton = 46;
    (list[131]).neutron = 56;
    (list[131]).mass_of_atom = 101.905608;
    (list[131]).mass_of_nucleus = 185716.47021271687;
    (list[131]).RAM_of_element = 106.42;
    (list[131]).abundency = 0.0102;
    STR_NAME_COPY((list[131]).symbol, "Pd");
    STR_NAME_COPY((list[131]).name, "Palladium");
    //===========================================================================
    (list[132]).proton = 46;
    (list[132]).neutron = 58;
    (list[132]).mass_of_atom = 103.904035;
    (list[132]).mass_of_nucleus = 189359.37803050634;
    (list[132]).RAM_of_element = 106.42;
    (list[132]).abundency = 0.1114;
    STR_NAME_COPY((list[132]).symbol, "Pd");
    STR_NAME_COPY((list[132]).name, "Palladium");
    //===========================================================================
    (list[133]).proton = 46;
    (list[133]).neutron = 59;
    (list[133]).mass_of_atom = 104.905084;
    (list[133]).mass_of_nucleus = 191184.17784960923;
    (list[133]).RAM_of_element = 106.42;
    (list[133]).abundency = 0.22329999999999997;
    STR_NAME_COPY((list[133]).symbol, "Pd");
    STR_NAME_COPY((list[133]).name, "Palladium");
    //===========================================================================
    (list[134]).proton = 46;
    (list[134]).neutron = 60;
    (list[134]).mass_of_atom = 105.903483;
    (list[134]).mass_of_nucleus = 193004.1470165456;
    (list[134]).RAM_of_element = 106.42;
    (list[134]).abundency = 0.2733;
    STR_NAME_COPY((list[134]).symbol, "Pd");
    STR_NAME_COPY((list[134]).name, "Palladium");
    //===========================================================================
    (list[135]).proton = 46;
    (list[135]).neutron = 62;
    (list[135]).mass_of_atom = 107.903894;
    (list[135]).mass_of_nucleus = 196650.67144335332;
    (list[135]).RAM_of_element = 106.42;
    (list[135]).abundency = 0.2646;
    STR_NAME_COPY((list[135]).symbol, "Pd");
    STR_NAME_COPY((list[135]).name, "Palladium");
    //===========================================================================
    (list[136]).proton = 46;
    (list[136]).neutron = 64;
    (list[136]).mass_of_atom = 109.905152;
    (list[136]).mass_of_nucleus = 200298.73985596673;
    (list[136]).RAM_of_element = 106.42;
    (list[136]).abundency = 0.11720000000000001;
    STR_NAME_COPY((list[136]).symbol, "Pd");
    STR_NAME_COPY((list[136]).name, "Palladium");
    //===========================================================================
    (list[137]).proton = 47;
    (list[137]).neutron = 60;
    (list[137]).mass_of_atom = 106.905093;
    (list[137]).mass_of_nucleus = 194828.96947559773;
    (list[137]).RAM_of_element = 107.8682;
    (list[137]).abundency = 0.51839;
    STR_NAME_COPY((list[137]).symbol, "Ag");
    STR_NAME_COPY((list[137]).name, "Silver");
    //===========================================================================
    (list[138]).proton = 47;
    (list[138]).neutron = 62;
    (list[138]).mass_of_atom = 108.904756;
    (list[138]).mass_of_nucleus = 198474.13038247317;
    (list[138]).RAM_of_element = 107.8682;
    (list[138]).abundency = 0.48161000000000004;
    STR_NAME_COPY((list[138]).symbol, "Ag");
    STR_NAME_COPY((list[138]).name, "Silver");
    //===========================================================================
    (list[139]).proton = 48;
    (list[139]).neutron = 58;
    (list[139]).mass_of_atom = 105.906458;
    (list[139]).mass_of_nucleus = 193007.57010718537;
    (list[139]).RAM_of_element = 112.411;
    (list[139]).abundency = 0.0125;
    STR_NAME_COPY((list[139]).symbol, "Cd");
    STR_NAME_COPY((list[139]).name, "Cadmium");
    //===========================================================================
    (list[140]).proton = 48;
    (list[140]).neutron = 60;
    (list[140]).mass_of_atom = 107.904183;
    (list[140]).mass_of_nucleus = 196649.19825787263;
    (list[140]).RAM_of_element = 112.411;
    (list[140]).abundency = 0.0089;
    STR_NAME_COPY((list[140]).symbol, "Cd");
    STR_NAME_COPY((list[140]).name, "Cadmium");
    //===========================================================================
    (list[141]).proton = 48;
    (list[141]).neutron = 62;
    (list[141]).mass_of_atom = 109.903006;
    (list[141]).mass_of_nucleus = 200292.82793915566;
    (list[141]).RAM_of_element = 112.411;
    (list[141]).abundency = 0.1249;
    STR_NAME_COPY((list[141]).symbol, "Cd");
    STR_NAME_COPY((list[141]).name, "Cadmium");
    //===========================================================================
    (list[142]).proton = 48;
    (list[142]).neutron = 63;
    (list[142]).mass_of_atom = 110.904182;
    (list[142]).mass_of_nucleus = 202117.85926498502;
    (list[142]).RAM_of_element = 112.411;
    (list[142]).abundency = 0.128;
    STR_NAME_COPY((list[142]).symbol, "Cd");
    STR_NAME_COPY((list[142]).name, "Cadmium");
    //===========================================================================
    (list[143]).proton = 48;
    (list[143]).neutron = 64;
    (list[143]).mass_of_atom = 111.902757;
    (list[143]).mass_of_nucleus = 203938.14926014075;
    (list[143]).RAM_of_element = 112.411;
    (list[143]).abundency = 0.2413;
    STR_NAME_COPY((list[143]).symbol, "Cd");
    STR_NAME_COPY((list[143]).name, "Cadmium");
    //===========================================================================
    (list[144]).proton = 48;
    (list[144]).neutron = 65;
    (list[144]).mass_of_atom = 112.904401;
    (list[144]).mass_of_nucleus = 205764.0336973716;
    (list[144]).RAM_of_element = 112.411;
    (list[144]).abundency = 0.1222;
    STR_NAME_COPY((list[144]).symbol, "Cd");
    STR_NAME_COPY((list[144]).name, "Cadmium");
    //===========================================================================
    (list[145]).proton = 48;
    (list[145]).neutron = 66;
    (list[145]).mass_of_atom = 113.903358;
    (list[145]).mass_of_nucleus = 207585.02003559438;
    (list[145]).RAM_of_element = 112.411;
    (list[145]).abundency = 0.2873;
    STR_NAME_COPY((list[145]).symbol, "Cd");
    STR_NAME_COPY((list[145]).name, "Cadmium");
    //===========================================================================
    (list[146]).proton = 48;
    (list[146]).neutron = 68;
    (list[146]).mass_of_atom = 115.904755;
    (list[146]).mass_of_nucleus = 211233.34182958555;
    (list[146]).RAM_of_element = 112.411;
    (list[146]).abundency = 0.07490000000000001;
    STR_NAME_COPY((list[146]).symbol, "Cd");
    STR_NAME_COPY((list[146]).name, "Cadmium");
    //===========================================================================
    (list[147]).proton = 49;
    (list[147]).neutron = 64;
    (list[147]).mass_of_atom = 112.904061;
    (list[147]).mass_of_nucleus = 205762.4139155842;
    (list[147]).RAM_of_element = 114.818;
    (list[147]).abundency = 0.0429;
    STR_NAME_COPY((list[147]).symbol, "In");
    STR_NAME_COPY((list[147]).name, "Indium");
    //===========================================================================
    (list[148]).proton = 49;
    (list[148]).neutron = 66;
    (list[148]).mass_of_atom = 114.903878;
    (list[148]).mass_of_nucleus = 209407.8555471516;
    (list[148]).RAM_of_element = 114.818;
    (list[148]).abundency = 0.9571;
    STR_NAME_COPY((list[148]).symbol, "In");
    STR_NAME_COPY((list[148]).name, "Indium");
    //===========================================================================
    (list[149]).proton = 50;
    (list[149]).neutron = 62;
    (list[149]).mass_of_atom = 111.904821;
    (list[149]).mass_of_nucleus = 203939.9117001678;
    (list[149]).RAM_of_element = 118.71;
    (list[149]).abundency = 0.0097;
    STR_NAME_COPY((list[149]).symbol, "Sn");
    STR_NAME_COPY((list[149]).name, "Tin");
    //===========================================================================
    (list[150]).proton = 50;
    (list[150]).neutron = 64;
    (list[150]).mass_of_atom = 113.902782;
    (list[150]).mass_of_nucleus = 207581.97005233102;
    (list[150]).RAM_of_element = 118.71;
    (list[150]).abundency = 0.0066;
    STR_NAME_COPY((list[150]).symbol, "Sn");
    STR_NAME_COPY((list[150]).name, "Tin");
    //===========================================================================
    (list[151]).proton = 50;
    (list[151]).neutron = 65;
    (list[151]).mass_of_atom = 114.903346;
    (list[151]).mass_of_nucleus = 209405.88577094304;
    (list[151]).RAM_of_element = 118.71;
    (list[151]).abundency = 0.0034000000000000002;
    STR_NAME_COPY((list[151]).symbol, "Sn");
    STR_NAME_COPY((list[151]).name, "Tin");
    //===========================================================================
    (list[152]).proton = 50;
    (list[152]).neutron = 66;
    (list[152]).mass_of_atom = 115.901744;
    (list[152]).mass_of_nucleus = 211225.85311499183;
    (list[152]).RAM_of_element = 118.71;
    (list[152]).abundency = 0.1454;
    STR_NAME_COPY((list[152]).symbol, "Sn");
    STR_NAME_COPY((list[152]).name, "Tin");
    //===========================================================================
    (list[153]).proton = 50;
    (list[153]).neutron = 67;
    (list[153]).mass_of_atom = 116.902954;
    (list[153]).mass_of_nucleus = 213050.94641899993;
    (list[153]).RAM_of_element = 118.71;
    (list[153]).abundency = 0.0768;
    STR_NAME_COPY((list[153]).symbol, "Sn");
    STR_NAME_COPY((list[153]).name, "Tin");
    //===========================================================================
    (list[154]).proton = 50;
    (list[154]).neutron = 68;
    (list[154]).mass_of_atom = 117.901606;
    (list[154]).mass_of_nucleus = 214871.37677650165;
    (list[154]).RAM_of_element = 118.71;
    (list[154]).abundency = 0.2422;
    STR_NAME_COPY((list[154]).symbol, "Sn");
    STR_NAME_COPY((list[154]).name, "Tin");
    //===========================================================================
    (list[155]).proton = 50;
    (list[155]).neutron = 69;
    (list[155]).mass_of_atom = 118.903309;
    (list[155]).mass_of_nucleus = 216697.36876410147;
    (list[155]).RAM_of_element = 118.71;
    (list[155]).abundency = 0.0859;
    STR_NAME_COPY((list[155]).symbol, "Sn");
    STR_NAME_COPY((list[155]).name, "Tin");
    //===========================================================================
    (list[156]).proton = 50;
    (list[156]).neutron = 70;
    (list[156]).mass_of_atom = 119.902197;
    (list[156]).mass_of_nucleus = 218518.22932307917;
    (list[156]).RAM_of_element = 118.71;
    (list[156]).abundency = 0.3258;
    STR_NAME_COPY((list[156]).symbol, "Sn");
    STR_NAME_COPY((list[156]).name, "Tin");
    //===========================================================================
    (list[157]).proton = 50;
    (list[157]).neutron = 72;
    (list[157]).mass_of_atom = 121.903440;
    (list[157]).mass_of_nucleus = 222166.2703923784;
    (list[157]).RAM_of_element = 118.71;
    (list[157]).abundency = 0.0463;
    STR_NAME_COPY((list[157]).symbol, "Sn");
    STR_NAME_COPY((list[157]).name, "Tin");
    //===========================================================================
    (list[158]).proton = 50;
    (list[158]).neutron = 74;
    (list[158]).mass_of_atom = 123.905275;
    (list[158]).mass_of_nucleus = 225815.39061114276;
    (list[158]).RAM_of_element = 118.71;
    (list[158]).abundency = 0.0579;
    STR_NAME_COPY((list[158]).symbol, "Sn");
    STR_NAME_COPY((list[158]).name, "Tin");
    //===========================================================================
    (list[159]).proton = 51;
    (list[159]).neutron = 70;
    (list[159]).mass_of_atom = 120.903818;
    (list[159]).mass_of_nucleus = 220343.071833895;
    (list[159]).RAM_of_element = 121.76;
    (list[159]).abundency = 0.5721;
    STR_NAME_COPY((list[159]).symbol, "Sb");
    STR_NAME_COPY((list[159]).name, "Antimony");
    //===========================================================================
    (list[160]).proton = 51;
    (list[160]).neutron = 72;
    (list[160]).mass_of_atom = 122.904216;
    (list[160]).mass_of_nucleus = 223989.57256316376;
    (list[160]).RAM_of_element = 121.76;
    (list[160]).abundency = 0.4279;
    STR_NAME_COPY((list[160]).symbol, "Sb");
    STR_NAME_COPY((list[160]).name, "Antimony");
    //===========================================================================
    (list[161]).proton = 52;
    (list[161]).neutron = 68;
    (list[161]).mass_of_atom = 119.904020;
    (list[161]).mass_of_nucleus = 218519.5524471922;
    (list[161]).RAM_of_element = 127.6;
    (list[161]).abundency = 0.0009;
    STR_NAME_COPY((list[161]).symbol, "Te");
    STR_NAME_COPY((list[161]).name, "Tellurium");
    //===========================================================================
    (list[162]).proton = 52;
    (list[162]).neutron = 70;
    (list[162]).mass_of_atom = 121.903047;
    (list[162]).mass_of_nucleus = 222163.55399754766;
    (list[162]).RAM_of_element = 127.6;
    (list[162]).abundency = 0.0255;
    STR_NAME_COPY((list[162]).symbol, "Te");
    STR_NAME_COPY((list[162]).name, "Tellurium");
    //===========================================================================
    (list[163]).proton = 52;
    (list[163]).neutron = 71;
    (list[163]).mass_of_atom = 122.904273;
    (list[163]).mass_of_nucleus = 223988.67646775753;
    (list[163]).RAM_of_element = 127.6;
    (list[163]).abundency = 0.0089;
    STR_NAME_COPY((list[163]).symbol, "Te");
    STR_NAME_COPY((list[163]).name, "Tellurium");
    //===========================================================================
    (list[164]).proton = 52;
    (list[164]).neutron = 72;
    (list[164]).mass_of_atom = 123.902819;
    (list[164]).mass_of_nucleus = 225808.91359917258;
    (list[164]).RAM_of_element = 127.6;
    (list[164]).abundency = 0.047400000000000005;
    STR_NAME_COPY((list[164]).symbol, "Te");
    STR_NAME_COPY((list[164]).name, "Tellurium");
    //===========================================================================
    (list[165]).proton = 52;
    (list[165]).neutron = 73;
    (list[165]).mass_of_atom = 124.904425;
    (list[165]).mass_of_nucleus = 227634.72876667426;
    (list[165]).RAM_of_element = 127.6;
    (list[165]).abundency = 0.0707;
    STR_NAME_COPY((list[165]).symbol, "Te");
    STR_NAME_COPY((list[165]).name, "Tellurium");
    //===========================================================================
    (list[166]).proton = 52;
    (list[166]).neutron = 74;
    (list[166]).mass_of_atom = 125.903306;
    (list[166]).mass_of_nucleus = 229455.57656543865;
    (list[166]).RAM_of_element = 127.6;
    (list[166]).abundency = 0.1884;
    STR_NAME_COPY((list[166]).symbol, "Te");
    STR_NAME_COPY((list[166]).name, "Tellurium");
    //===========================================================================
    (list[167]).proton = 52;
    (list[167]).neutron = 76;
    (list[167]).mass_of_atom = 127.904461;
    (list[167]).mass_of_nucleus = 233103.45722062822;
    (list[167]).RAM_of_element = 127.6;
    (list[167]).abundency = 0.31739999999999996;
    STR_NAME_COPY((list[167]).symbol, "Te");
    STR_NAME_COPY((list[167]).name, "Tellurium");
    //===========================================================================
    (list[168]).proton = 52;
    (list[168]).neutron = 78;
    (list[168]).mass_of_atom = 129.906223;
    (list[168]).mass_of_nucleus = 236752.44436859706;
    (list[168]).RAM_of_element = 127.6;
    (list[168]).abundency = 0.3408;
    STR_NAME_COPY((list[168]).symbol, "Te");
    STR_NAME_COPY((list[168]).name, "Tellurium");
    //===========================================================================
    (list[169]).proton = 53;
    (list[169]).neutron = 74;
    (list[169]).mass_of_atom = 126.904468;
    (list[169]).mass_of_nucleus = 231279.58237084147;
    (list[169]).RAM_of_element = 126.9045;
    (list[169]).abundency = 1.0;
    STR_NAME_COPY((list[169]).symbol, "I");
    STR_NAME_COPY((list[169]).name, "Iodine");
    //===========================================================================
    (list[170]).proton = 54;
    (list[170]).neutron = 70;
    (list[170]).mass_of_atom = 123.905896;
    (list[170]).mass_of_nucleus = 225812.52262434855;
    (list[170]).RAM_of_element = 131.293;
    (list[170]).abundency = 0.0009;
    STR_NAME_COPY((list[170]).symbol, "Xe");
    STR_NAME_COPY((list[170]).name, "Xenon");
    //===========================================================================
    (list[171]).proton = 54;
    (list[171]).neutron = 72;
    (list[171]).mass_of_atom = 125.904269;
    (list[171]).mass_of_nucleus = 229455.33200620709;
    (list[171]).RAM_of_element = 131.293;
    (list[171]).abundency = 0.0009;
    STR_NAME_COPY((list[171]).symbol, "Xe");
    STR_NAME_COPY((list[171]).name, "Xenon");
    //===========================================================================
    (list[172]).proton = 54;
    (list[172]).neutron = 74;
    (list[172]).mass_of_atom = 127.903530;
    (list[172]).mass_of_nucleus = 233099.7601122633;
    (list[172]).RAM_of_element = 131.293;
    (list[172]).abundency = 0.0192;
    STR_NAME_COPY((list[172]).symbol, "Xe");
    STR_NAME_COPY((list[172]).name, "Xenon");
    //===========================================================================
    (list[173]).proton = 54;
    (list[173]).neutron = 75;
    (list[173]).mass_of_atom = 128.904779;
    (list[173]).mass_of_nucleus = 234924.92450888816;
    (list[173]).RAM_of_element = 131.293;
    (list[173]).abundency = 0.2644;
    STR_NAME_COPY((list[173]).symbol, "Xe");
    STR_NAME_COPY((list[173]).name, "Xenon");
    //===========================================================================
    (list[174]).proton = 54;
    (list[174]).neutron = 76;
    (list[174]).mass_of_atom = 129.903508;
    (list[174]).mass_of_nucleus = 236745.49522873585;
    (list[174]).RAM_of_element = 131.293;
    (list[174]).abundency = 0.0408;
    STR_NAME_COPY((list[174]).symbol, "Xe");
    STR_NAME_COPY((list[174]).name, "Xenon");
    //===========================================================================
    (list[175]).proton = 54;
    (list[175]).neutron = 77;
    (list[175]).mass_of_atom = 130.905082;
    (list[175]).mass_of_nucleus = 238571.25206383402;
    (list[175]).RAM_of_element = 131.293;
    (list[175]).abundency = 0.2118;
    STR_NAME_COPY((list[175]).symbol, "Xe");
    STR_NAME_COPY((list[175]).name, "Xenon");
    //===========================================================================
    (list[176]).proton = 54;
    (list[176]).neutron = 78;
    (list[176]).mass_of_atom = 131.904154;
    (list[176]).mass_of_nucleus = 240392.44803413196;
    (list[176]).RAM_of_element = 131.293;
    (list[176]).abundency = 0.26890000000000003;
    STR_NAME_COPY((list[176]).symbol, "Xe");
    STR_NAME_COPY((list[176]).name, "Xenon");
    //===========================================================================
    (list[177]).proton = 54;
    (list[177]).neutron = 80;
    (list[177]).mass_of_atom = 133.905395;
    (list[177]).mass_of_nucleus = 244040.48545765594;
    (list[177]).RAM_of_element = 131.293;
    (list[177]).abundency = 0.10439999999999999;
    STR_NAME_COPY((list[177]).symbol, "Xe");
    STR_NAME_COPY((list[177]).name, "Xenon");
    //===========================================================================
    (list[178]).proton = 54;
    (list[178]).neutron = 82;
    (list[178]).mass_of_atom = 135.907220;
    (list[178]).mass_of_nucleus = 247689.58744754418;
    (list[178]).RAM_of_element = 131.293;
    (list[178]).abundency = 0.08869999999999999;
    STR_NAME_COPY((list[178]).symbol, "Xe");
    STR_NAME_COPY((list[178]).name, "Xenon");
    //===========================================================================
    (list[179]).proton = 55;
    (list[179]).neutron = 78;
    (list[179]).mass_of_atom = 132.905447;
    (list[179]).mass_of_nucleus = 242216.69263781168;
    (list[179]).RAM_of_element = 132.9055;
    (list[179]).abundency = 1.0;
    STR_NAME_COPY((list[179]).symbol, "Cs");
    STR_NAME_COPY((list[179]).name, "Cesium");
    //===========================================================================
    (list[180]).proton = 56;
    (list[180]).neutron = 74;
    (list[180]).mass_of_atom = 129.906310;
    (list[180]).mass_of_nucleus = 236748.60295981908;
    (list[180]).RAM_of_element = 137.327;
    (list[180]).abundency = 0.00106;
    STR_NAME_COPY((list[180]).symbol, "Ba");
    STR_NAME_COPY((list[180]).name, "Barium");
    //===========================================================================
    (list[181]).proton = 56;
    (list[181]).neutron = 76;
    (list[181]).mass_of_atom = 131.905056;
    (list[181]).mass_of_nucleus = 240392.09227875617;
    (list[181]).RAM_of_element = 137.327;
    (list[181]).abundency = 0.00101;
    STR_NAME_COPY((list[181]).symbol, "Ba");
    STR_NAME_COPY((list[181]).name, "Barium");
    //===========================================================================
    (list[182]).proton = 56;
    (list[182]).neutron = 78;
    (list[182]).mass_of_atom = 133.904503;
    (list[182]).mass_of_nucleus = 244036.85944190784;
    (list[182]).RAM_of_element = 137.327;
    (list[182]).abundency = 0.024169999999999997;
    STR_NAME_COPY((list[182]).symbol, "Ba");
    STR_NAME_COPY((list[182]).name, "Barium");
    //===========================================================================
    (list[183]).proton = 56;
    (list[183]).neutron = 79;
    (list[183]).mass_of_atom = 134.905683;
    (list[183]).mass_of_nucleus = 245861.89805928766;
    (list[183]).RAM_of_element = 137.327;
    (list[183]).abundency = 0.06591999999999999;
    STR_NAME_COPY((list[183]).symbol, "Ba");
    STR_NAME_COPY((list[183]).name, "Barium");
    //===========================================================================
    (list[184]).proton = 56;
    (list[184]).neutron = 80;
    (list[184]).mass_of_atom = 135.904570;
    (list[184]).mass_of_nucleus = 247682.7567953777;
    (list[184]).RAM_of_element = 137.327;
    (list[184]).abundency = 0.07854;
    STR_NAME_COPY((list[184]).symbol, "Ba");
    STR_NAME_COPY((list[184]).name, "Barium");
    //===========================================================================
    (list[185]).proton = 56;
    (list[185]).neutron = 81;
    (list[185]).mass_of_atom = 136.905821;
    (list[185]).mass_of_nucleus = 249507.92483777783;
    (list[185]).RAM_of_element = 137.327;
    (list[185]).abundency = 0.11231999999999999;
    STR_NAME_COPY((list[185]).symbol, "Ba");
    STR_NAME_COPY((list[185]).name, "Barium");
    //===========================================================================
    (list[186]).proton = 56;
    (list[186]).neutron = 82;
    (list[186]).mass_of_atom = 137.905241;
    (list[186]).mass_of_nucleus = 251329.75517296398;
    (list[186]).RAM_of_element = 137.327;
    (list[186]).abundency = 0.71698;
    STR_NAME_COPY((list[186]).symbol, "Ba");
    STR_NAME_COPY((list[186]).name, "Barium");
    //===========================================================================
    (list[187]).proton = 57;
    (list[187]).neutron = 81;
    (list[187]).mass_of_atom = 137.907107;
    (list[187]).mass_of_nucleus = 251332.15668124426;
    (list[187]).RAM_of_element = 138.9055;
    (list[187]).abundency = 0.0009;
    STR_NAME_COPY((list[187]).symbol, "La");
    STR_NAME_COPY((list[187]).name, "Lanthanum");
    //===========================================================================
    (list[188]).proton = 57;
    (list[188]).neutron = 82;
    (list[188]).mass_of_atom = 138.906348;
    (list[188]).mass_of_nucleus = 253153.6607195483;
    (list[188]).RAM_of_element = 138.9055;
    (list[188]).abundency = 0.9991;
    STR_NAME_COPY((list[188]).symbol, "La");
    STR_NAME_COPY((list[188]).name, "Lanthanum");
    //===========================================================================
    (list[189]).proton = 58;
    (list[189]).neutron = 78;
    (list[189]).mass_of_atom = 135.907144;
    (list[189]).mass_of_nucleus = 247685.44890808582;
    (list[189]).RAM_of_element = 140.116;
    (list[189]).abundency = 0.00185;
    STR_NAME_COPY((list[189]).symbol, "Ce");
    STR_NAME_COPY((list[189]).name, "Cerium");
    //===========================================================================
    (list[190]).proton = 58;
    (list[190]).neutron = 80;
    (list[190]).mass_of_atom = 137.905986;
    (list[190]).mass_of_nucleus = 251329.11322423347;
    (list[190]).RAM_of_element = 140.116;
    (list[190]).abundency = 0.00251;
    STR_NAME_COPY((list[190]).symbol, "Ce");
    STR_NAME_COPY((list[190]).name, "Cerium");
    //===========================================================================
    (list[191]).proton = 58;
    (list[191]).neutron = 82;
    (list[191]).mass_of_atom = 139.905434;
    (list[191]).mass_of_nucleus = 254973.88221027277;
    (list[191]).RAM_of_element = 140.116;
    (list[191]).abundency = 0.8845000000000001;
    STR_NAME_COPY((list[191]).symbol, "Ce");
    STR_NAME_COPY((list[191]).name, "Cerium");
    //===========================================================================
    (list[192]).proton = 58;
    (list[192]).neutron = 84;
    (list[192]).mass_of_atom = 141.909240;
    (list[192]).mass_of_nucleus = 258626.59534051642;
    (list[192]).RAM_of_element = 140.116;
    (list[192]).abundency = 0.11114;
    STR_NAME_COPY((list[192]).symbol, "Ce");
    STR_NAME_COPY((list[192]).name, "Cerium");
    //===========================================================================
    (list[193]).proton = 59;
    (list[193]).neutron = 82;
    (list[193]).mass_of_atom = 140.907648;
    (list[193]).mass_of_nucleus = 256799.80569344127;
    (list[193]).RAM_of_element = 140.9077;
    (list[193]).abundency = 1.0;
    STR_NAME_COPY((list[193]).symbol, "Pr");
    STR_NAME_COPY((list[193]).name, "Praseodymium");
    //===========================================================================
    (list[194]).proton = 60;
    (list[194]).neutron = 82;
    (list[194]).mass_of_atom = 141.907719;
    (list[194]).mass_of_nucleus = 258621.82272846156;
    (list[194]).RAM_of_element = 144.24;
    (list[194]).abundency = 0.272;
    STR_NAME_COPY((list[194]).symbol, "Nd");
    STR_NAME_COPY((list[194]).name, "Neodymium");
    //===========================================================================
    (list[195]).proton = 60;
    (list[195]).neutron = 83;
    (list[195]).mass_of_atom = 142.909810;
    (list[195]).mass_of_nucleus = 260448.52199645407;
    (list[195]).RAM_of_element = 144.24;
    (list[195]).abundency = 0.122;
    STR_NAME_COPY((list[195]).symbol, "Nd");
    STR_NAME_COPY((list[195]).name, "Neodymium");
    //===========================================================================
    (list[196]).proton = 60;
    (list[196]).neutron = 84;
    (list[196]).mass_of_atom = 143.910083;
    (list[196]).mass_of_nucleus = 262271.9072547716;
    (list[196]).RAM_of_element = 144.24;
    (list[196]).abundency = 0.23800000000000002;
    STR_NAME_COPY((list[196]).symbol, "Nd");
    STR_NAME_COPY((list[196]).name, "Neodymium");
    //===========================================================================
    (list[197]).proton = 60;
    (list[197]).neutron = 85;
    (list[197]).mass_of_atom = 144.912569;
    (list[197]).mass_of_nucleus = 264099.32656337006;
    (list[197]).RAM_of_element = 144.24;
    (list[197]).abundency = 0.083;
    STR_NAME_COPY((list[197]).symbol, "Nd");
    STR_NAME_COPY((list[197]).name, "Neodymium");
    //===========================================================================
    (list[198]).proton = 60;
    (list[198]).neutron = 86;
    (list[198]).mass_of_atom = 145.913112;
    (list[198]).mass_of_nucleus = 265923.20400134235;
    (list[198]).RAM_of_element = 144.24;
    (list[198]).abundency = 0.172;
    STR_NAME_COPY((list[198]).symbol, "Nd");
    STR_NAME_COPY((list[198]).name, "Neodymium");
    //===========================================================================
    (list[199]).proton = 60;
    (list[199]).neutron = 88;
    (list[199]).mass_of_atom = 147.916889;
    (list[199]).mass_of_nucleus = 269575.8642678453;
    (list[199]).RAM_of_element = 144.24;
    (list[199]).abundency = 0.057;
    STR_NAME_COPY((list[199]).symbol, "Nd");
    STR_NAME_COPY((list[199]).name, "Neodymium");
    //===========================================================================
    (list[200]).proton = 60;
    (list[200]).neutron = 90;
    (list[200]).mass_of_atom = 149.920887;
    (list[200]).mass_of_nucleus = 273228.92739251006;
    (list[200]).RAM_of_element = 144.24;
    (list[200]).abundency = 0.055999999999999994;
    STR_NAME_COPY((list[200]).symbol, "Nd");
    STR_NAME_COPY((list[200]).name, "Neodymium");
    //===========================================================================
    (list[201]).proton = 61;
    (list[201]).neutron = 84;
    (list[201]).mass_of_atom = 144.912744;
    (list[201]).mass_of_nucleus = 264098.64556870185;
    (list[201]).RAM_of_element = 145;
    (list[201]).abundency = 0.0;
    STR_NAME_COPY((list[201]).symbol, "Pm");
    STR_NAME_COPY((list[201]).name, "Promethium");
    //===========================================================================
    (list[202]).proton = 62;
    (list[202]).neutron = 82;
    (list[202]).mass_of_atom = 143.911995;
    (list[202]).mass_of_nucleus = 262273.3926158819;
    (list[202]).RAM_of_element = 150.36;
    (list[202]).abundency = 0.030699999999999998;
    STR_NAME_COPY((list[202]).symbol, "Sm");
    STR_NAME_COPY((list[202]).name, "Samarium");
    //===========================================================================
    (list[203]).proton = 62;
    (list[203]).neutron = 85;
    (list[203]).mass_of_atom = 146.914893;
    (list[203]).mass_of_nucleus = 267747.33817417576;
    (list[203]).RAM_of_element = 150.36;
    (list[203]).abundency = 0.1499;
    STR_NAME_COPY((list[203]).symbol, "Sm");
    STR_NAME_COPY((list[203]).name, "Samarium");
    //===========================================================================
    (list[204]).proton = 62;
    (list[204]).neutron = 86;
    (list[204]).mass_of_atom = 147.914818;
    (list[204]).mass_of_nucleus = 269570.089067605;
    (list[204]).RAM_of_element = 150.36;
    (list[204]).abundency = 0.1124;
    STR_NAME_COPY((list[204]).symbol, "Sm");
    STR_NAME_COPY((list[204]).name, "Samarium");
    //===========================================================================
    (list[205]).proton = 62;
    (list[205]).neutron = 87;
    (list[205]).mass_of_atom = 148.917180;
    (list[205]).mass_of_nucleus = 271397.2823381398;
    (list[205]).RAM_of_element = 150.36;
    (list[205]).abundency = 0.1382;
    STR_NAME_COPY((list[205]).symbol, "Sm");
    STR_NAME_COPY((list[205]).name, "Samarium");
    //===========================================================================
    (list[206]).proton = 62;
    (list[206]).neutron = 88;
    (list[206]).mass_of_atom = 149.917271;
    (list[206]).mass_of_nucleus = 273220.3358309123;
    (list[206]).RAM_of_element = 150.36;
    (list[206]).abundency = 0.0738;
    STR_NAME_COPY((list[206]).symbol, "Sm");
    STR_NAME_COPY((list[206]).name, "Samarium");
    //===========================================================================
    (list[207]).proton = 62;
    (list[207]).neutron = 90;
    (list[207]).mass_of_atom = 151.919728;
    (list[207]).mass_of_nucleus = 276870.58988577005;
    (list[207]).RAM_of_element = 150.36;
    (list[207]).abundency = 0.2675;
    STR_NAME_COPY((list[207]).symbol, "Sm");
    STR_NAME_COPY((list[207]).name, "Samarium");
    //===========================================================================
    (list[208]).proton = 62;
    (list[208]).neutron = 92;
    (list[208]).mass_of_atom = 153.922205;
    (list[208]).mass_of_nucleus = 280520.88039838005;
    (list[208]).RAM_of_element = 150.36;
    (list[208]).abundency = 0.2275;
    STR_NAME_COPY((list[208]).symbol, "Sm");
    STR_NAME_COPY((list[208]).name, "Samarium");
    //===========================================================================
    (list[209]).proton = 63;
    (list[209]).neutron = 88;
    (list[209]).mass_of_atom = 150.919846;
    (list[209]).mass_of_nucleus = 275046.91737650806;
    (list[209]).RAM_of_element = 151.964;
    (list[209]).abundency = 0.4781;
    STR_NAME_COPY((list[209]).symbol, "Eu");
    STR_NAME_COPY((list[209]).name, "Europium");
    //===========================================================================
    (list[210]).proton = 63;
    (list[210]).neutron = 90;
    (list[210]).mass_of_atom = 152.921226;
    (list[210]).mass_of_nucleus = 278695.20818140986;
    (list[210]).RAM_of_element = 151.964;
    (list[210]).abundency = 0.5219;
    STR_NAME_COPY((list[210]).symbol, "Eu");
    STR_NAME_COPY((list[210]).name, "Europium");
    //===========================================================================
    (list[211]).proton = 64;
    (list[211]).neutron = 88;
    (list[211]).mass_of_atom = 151.919788;
    (list[211]).mass_of_nucleus = 276868.6992590267;
    (list[211]).RAM_of_element = 157.25;
    (list[211]).abundency = 0.002;
    STR_NAME_COPY((list[211]).symbol, "Gd");
    STR_NAME_COPY((list[211]).name, "Gadolinium");
    //===========================================================================
    (list[212]).proton = 64;
    (list[212]).neutron = 90;
    (list[212]).mass_of_atom = 153.920862;
    (list[212]).mass_of_nucleus = 280516.4322603198;
    (list[212]).RAM_of_element = 157.25;
    (list[212]).abundency = 0.0218;
    STR_NAME_COPY((list[212]).symbol, "Gd");
    STR_NAME_COPY((list[212]).name, "Gadolinium");
    //===========================================================================
    (list[213]).proton = 64;
    (list[213]).neutron = 91;
    (list[213]).mass_of_atom = 154.922619;
    (list[213]).mass_of_nucleus = 282342.5226838506;
    (list[213]).RAM_of_element = 157.25;
    (list[213]).abundency = 0.14800000000000002;
    STR_NAME_COPY((list[213]).symbol, "Gd");
    STR_NAME_COPY((list[213]).name, "Gadolinium");
    //===========================================================================
    (list[214]).proton = 64;
    (list[214]).neutron = 92;
    (list[214]).mass_of_atom = 155.922120;
    (list[214]).mass_of_nucleus = 284164.50067293324;
    (list[214]).RAM_of_element = 157.25;
    (list[214]).abundency = 0.2047;
    STR_NAME_COPY((list[214]).symbol, "Gd");
    STR_NAME_COPY((list[214]).name, "Gadolinium");
    //===========================================================================
    (list[215]).proton = 64;
    (list[215]).neutron = 93;
    (list[215]).mass_of_atom = 156.923957;
    (list[215]).mass_of_nucleus = 285990.73692747275;
    (list[215]).RAM_of_element = 157.25;
    (list[215]).abundency = 0.1565;
    STR_NAME_COPY((list[215]).symbol, "Gd");
    STR_NAME_COPY((list[215]).name, "Gadolinium");
    //===========================================================================
    (list[216]).proton = 64;
    (list[216]).neutron = 94;
    (list[216]).mass_of_atom = 157.924101;
    (list[216]).mass_of_nucleus = 287813.8870332886;
    (list[216]).RAM_of_element = 157.25;
    (list[216]).abundency = 0.2484;
    STR_NAME_COPY((list[216]).symbol, "Gd");
    STR_NAME_COPY((list[216]).name, "Gadolinium");
    //===========================================================================
    (list[217]).proton = 64;
    (list[217]).neutron = 96;
    (list[217]).mass_of_atom = 159.927051;
    (list[217]).mass_of_nucleus = 291465.0397717381;
    (list[217]).RAM_of_element = 157.25;
    (list[217]).abundency = 0.2186;
    STR_NAME_COPY((list[217]).symbol, "Gd");
    STR_NAME_COPY((list[217]).name, "Gadolinium");
    //===========================================================================
    (list[218]).proton = 65;
    (list[218]).neutron = 94;
    (list[218]).mass_of_atom = 158.925343;
    (list[218]).mass_of_nucleus = 289638.0386697002;
    (list[218]).RAM_of_element = 158.9253;
    (list[218]).abundency = 1.0;
    STR_NAME_COPY((list[218]).symbol, "Tb");
    STR_NAME_COPY((list[218]).name, "Terbium");
    //===========================================================================
    (list[219]).proton = 66;
    (list[219]).neutron = 90;
    (list[219]).mass_of_atom = 155.924278;
    (list[219]).mass_of_nucleus = 284166.43446439557;
    (list[219]).RAM_of_element = 162.5;
    (list[219]).abundency = 0.0006;
    STR_NAME_COPY((list[219]).symbol, "Dy");
    STR_NAME_COPY((list[219]).name, "Dysprosium");
    //===========================================================================
    (list[220]).proton = 66;
    (list[220]).neutron = 92;
    (list[220]).mass_of_atom = 157.924405;
    (list[220]).mass_of_nucleus = 287812.44119112205;
    (list[220]).RAM_of_element = 162.5;
    (list[220]).abundency = 0.001;
    STR_NAME_COPY((list[220]).symbol, "Dy");
    STR_NAME_COPY((list[220]).name, "Dysprosium");
    //===========================================================================
    (list[221]).proton = 66;
    (list[221]).neutron = 94;
    (list[221]).mass_of_atom = 159.925194;
    (list[221]).mass_of_nucleus = 291459.65466944635;
    (list[221]).RAM_of_element = 162.5;
    (list[221]).abundency = 0.023399999999999997;
    STR_NAME_COPY((list[221]).symbol, "Dy");
    STR_NAME_COPY((list[221]).name, "Dysprosium");
    //===========================================================================
    (list[222]).proton = 66;
    (list[222]).neutron = 95;
    (list[222]).mass_of_atom = 160.926930;
    (list[222]).mass_of_nucleus = 293285.7068123373;
    (list[222]).RAM_of_element = 162.5;
    (list[222]).abundency = 0.1891;
    STR_NAME_COPY((list[222]).symbol, "Dy");
    STR_NAME_COPY((list[222]).name, "Dysprosium");
    //===========================================================================
    (list[223]).proton = 66;
    (list[223]).neutron = 96;
    (list[223]).mass_of_atom = 161.926795;
    (list[223]).mass_of_nucleus = 295108.3483325099;
    (list[223]).RAM_of_element = 162.5;
    (list[223]).abundency = 0.2551;
    STR_NAME_COPY((list[223]).symbol, "Dy");
    STR_NAME_COPY((list[223]).name, "Dysprosium");
    //===========================================================================
    (list[224]).proton = 66;
    (list[224]).neutron = 97;
    (list[224]).mass_of_atom = 162.928728;
    (list[224]).mass_of_nucleus = 296934.7595842601;
    (list[224]).RAM_of_element = 162.5;
    (list[224]).abundency = 0.249;
    STR_NAME_COPY((list[224]).symbol, "Dy");
    STR_NAME_COPY((list[224]).name, "Dysprosium");
    //===========================================================================
    (list[225]).proton = 66;
    (list[225]).neutron = 98;
    (list[225]).mass_of_atom = 163.929171;
    (list[225]).mass_of_nucleus = 298758.4547334713;
    (list[225]).RAM_of_element = 162.5;
    (list[225]).abundency = 0.2818;
    STR_NAME_COPY((list[225]).symbol, "Dy");
    STR_NAME_COPY((list[225]).name, "Dysprosium");
    //===========================================================================
    (list[226]).proton = 67;
    (list[226]).neutron = 98;
    (list[226]).mass_of_atom = 164.930319;
    (list[226]).mass_of_nucleus = 300582.4350184476;
    (list[226]).RAM_of_element = 164.9303;
    (list[226]).abundency = 1.0;
    STR_NAME_COPY((list[226]).symbol, "Ho");
    STR_NAME_COPY((list[226]).name, "Holmium");
    //===========================================================================
    (list[227]).proton = 68;
    (list[227]).neutron = 94;
    (list[227]).mass_of_atom = 161.928775;
    (list[227]).mass_of_nucleus = 295109.95764997776;
    (list[227]).RAM_of_element = 167.259;
    (list[227]).abundency = 0.0014000000000000002;
    STR_NAME_COPY((list[227]).symbol, "Er");
    STR_NAME_COPY((list[227]).name, "Erbium");
    //===========================================================================
    (list[228]).proton = 68;
    (list[228]).neutron = 96;
    (list[228]).mass_of_atom = 163.929197;
    (list[228]).mass_of_nucleus = 298756.50212854915;
    (list[228]).RAM_of_element = 167.259;
    (list[228]).abundency = 0.0161;
    STR_NAME_COPY((list[228]).symbol, "Er");
    STR_NAME_COPY((list[228]).name, "Erbium");
    //===========================================================================
    (list[229]).proton = 68;
    (list[229]).neutron = 98;
    (list[229]).mass_of_atom = 165.930290;
    (list[229]).mass_of_nucleus = 302404.2697647069;
    (list[229]).RAM_of_element = 167.259;
    (list[229]).abundency = 0.3361;
    STR_NAME_COPY((list[229]).symbol, "Er");
    STR_NAME_COPY((list[229]).name, "Erbium");
    //===========================================================================
    (list[230]).proton = 68;
    (list[230]).neutron = 99;
    (list[230]).mass_of_atom = 166.932045;
    (list[230]).mass_of_nucleus = 304230.3565424624;
    (list[230]).RAM_of_element = 167.259;
    (list[230]).abundency = 0.2293;
    STR_NAME_COPY((list[230]).symbol, "Er");
    STR_NAME_COPY((list[230]).name, "Erbium");
    //===========================================================================
    (list[231]).proton = 68;
    (list[231]).neutron = 100;
    (list[231]).mass_of_atom = 167.932368;
    (list[231]).mass_of_nucleus = 306053.83294516045;
    (list[231]).RAM_of_element = 167.259;
    (list[231]).abundency = 0.26780000000000004;
    STR_NAME_COPY((list[231]).symbol, "Er");
    STR_NAME_COPY((list[231]).name, "Erbium");
    //===========================================================================
    (list[232]).proton = 68;
    (list[232]).neutron = 102;
    (list[232]).mass_of_atom = 169.935460;
    (list[232]).mass_of_nucleus = 309705.24453365064;
    (list[232]).RAM_of_element = 167.259;
    (list[232]).abundency = 0.1493;
    STR_NAME_COPY((list[232]).symbol, "Er");
    STR_NAME_COPY((list[232]).name, "Erbium");
    //===========================================================================
    (list[233]).proton = 69;
    (list[233]).neutron = 100;
    (list[233]).mass_of_atom = 168.934211;
    (list[233]).mass_of_nucleus = 307879.08013702574;
    (list[233]).RAM_of_element = 168.9342;
    (list[233]).abundency = 1.0;
    STR_NAME_COPY((list[233]).symbol, "Tm");
    STR_NAME_COPY((list[233]).name, "Thulium");
    //===========================================================================
    (list[234]).proton = 70;
    (list[234]).neutron = 98;
    (list[234]).mass_of_atom = 167.933894;
    (list[234]).mass_of_nucleus = 306054.61467165337;
    (list[234]).RAM_of_element = 173.04;
    (list[234]).abundency = 0.0013;
    STR_NAME_COPY((list[234]).symbol, "Yb");
    STR_NAME_COPY((list[234]).name, "Ytterbium");
    //===========================================================================
    (list[235]).proton = 70;
    (list[235]).neutron = 100;
    (list[235]).mass_of_atom = 169.934759;
    (list[235]).mass_of_nucleus = 309701.966689436;
    (list[235]).RAM_of_element = 173.04;
    (list[235]).abundency = 0.0304;
    STR_NAME_COPY((list[235]).symbol, "Yb");
    STR_NAME_COPY((list[235]).name, "Ytterbium");
    //===========================================================================
    (list[236]).proton = 70;
    (list[236]).neutron = 101;
    (list[236]).mass_of_atom = 170.936322;
    (list[236]).mass_of_nucleus = 311527.7034727704;
    (list[236]).RAM_of_element = 173.04;
    (list[236]).abundency = 0.14279999999999998;
    STR_NAME_COPY((list[236]).symbol, "Yb");
    STR_NAME_COPY((list[236]).name, "Ytterbium");
    //===========================================================================
    (list[237]).proton = 70;
    (list[237]).neutron = 102;
    (list[237]).mass_of_atom = 171.936378;
    (list[237]).mass_of_nucleus = 313350.69316447654;
    (list[237]).RAM_of_element = 173.04;
    (list[237]).abundency = 0.2183;
    STR_NAME_COPY((list[237]).symbol, "Yb");
    STR_NAME_COPY((list[237]).name, "Ytterbium");
    //===========================================================================
    (list[238]).proton = 70;
    (list[238]).neutron = 103;
    (list[238]).mass_of_atom = 172.938207;
    (list[238]).mass_of_nucleus = 315176.91483591526;
    (list[238]).RAM_of_element = 173.04;
    (list[238]).abundency = 0.1613;
    STR_NAME_COPY((list[238]).symbol, "Yb");
    STR_NAME_COPY((list[238]).name, "Ytterbium");
    //===========================================================================
    (list[239]).proton = 70;
    (list[239]).neutron = 104;
    (list[239]).mass_of_atom = 173.938858;
    (list[239]).mass_of_nucleus = 317000.9891457494;
    (list[239]).RAM_of_element = 173.04;
    (list[239]).abundency = 0.31829999999999997;
    STR_NAME_COPY((list[239]).symbol, "Yb");
    STR_NAME_COPY((list[239]).name, "Ytterbium");
    //===========================================================================
    (list[240]).proton = 70;
    (list[240]).neutron = 106;
    (list[240]).mass_of_atom = 175.942568;
    (list[240]).mass_of_nucleus = 320653.52727878245;
    (list[240]).RAM_of_element = 173.04;
    (list[240]).abundency = 0.1276;
    STR_NAME_COPY((list[240]).symbol, "Yb");
    STR_NAME_COPY((list[240]).name, "Ytterbium");
    //===========================================================================
    (list[241]).proton = 71;
    (list[241]).neutron = 104;
    (list[241]).mass_of_atom = 174.940768;
    (list[241]).mass_of_nucleus = 318826.35847108444;
    (list[241]).RAM_of_element = 174.967;
    (list[241]).abundency = 0.9741;
    STR_NAME_COPY((list[241]).symbol, "Lu");
    STR_NAME_COPY((list[241]).name, "Lutetium");
    //===========================================================================
    (list[242]).proton = 71;
    (list[242]).neutron = 105;
    (list[242]).mass_of_atom = 175.942682;
    (list[242]).mass_of_nucleus = 320652.73508797;
    (list[242]).RAM_of_element = 174.967;
    (list[242]).abundency = 0.0259;
    STR_NAME_COPY((list[242]).symbol, "Lu");
    STR_NAME_COPY((list[242]).name, "Lutetium");
    //===========================================================================
    (list[243]).proton = 72;
    (list[243]).neutron = 102;
    (list[243]).mass_of_atom = 173.940040;
    (list[243]).mass_of_nucleus = 317001.14379890444;
    (list[243]).RAM_of_element = 178.49;
    (list[243]).abundency = 0.0016;
    STR_NAME_COPY((list[243]).symbol, "Hf");
    STR_NAME_COPY((list[243]).name, "Hafnium");
    //===========================================================================
    (list[244]).proton = 72;
    (list[244]).neutron = 104;
    (list[244]).mass_of_atom = 175.941402;
    (list[244]).mass_of_nucleus = 320649.40179182927;
    (list[244]).RAM_of_element = 178.49;
    (list[244]).abundency = 0.0526;
    STR_NAME_COPY((list[244]).symbol, "Hf");
    STR_NAME_COPY((list[244]).name, "Hafnium");
    //===========================================================================
    (list[245]).proton = 72;
    (list[245]).neutron = 105;
    (list[245]).mass_of_atom = 176.943220;
    (list[245]).mass_of_nucleus = 322475.6034115042;
    (list[245]).RAM_of_element = 178.49;
    (list[245]).abundency = 0.18600000000000003;
    STR_NAME_COPY((list[245]).symbol, "Hf");
    STR_NAME_COPY((list[245]).name, "Hafnium");
    //===========================================================================
    (list[246]).proton = 72;
    (list[246]).neutron = 106;
    (list[246]).mass_of_atom = 177.943698;
    (list[246]).mass_of_nucleus = 324299.3623617818;
    (list[246]).RAM_of_element = 178.49;
    (list[246]).abundency = 0.2728;
    STR_NAME_COPY((list[246]).symbol, "Hf");
    STR_NAME_COPY((list[246]).name, "Hafnium");
    //===========================================================================
    (list[247]).proton = 72;
    (list[247]).neutron = 107;
    (list[247]).mass_of_atom = 178.945815;
    (list[247]).mass_of_nucleus = 326126.10902485217;
    (list[247]).RAM_of_element = 178.49;
    (list[247]).abundency = 0.1362;
    STR_NAME_COPY((list[247]).symbol, "Hf");
    STR_NAME_COPY((list[247]).name, "Hafnium");
    //===========================================================================
    (list[248]).proton = 72;
    (list[248]).neutron = 108;
    (list[248]).mass_of_atom = 179.946549;
    (list[248]).mass_of_nucleus = 327950.3346343579;
    (list[248]).RAM_of_element = 178.49;
    (list[248]).abundency = 0.3508;
    STR_NAME_COPY((list[248]).symbol, "Hf");
    STR_NAME_COPY((list[248]).name, "Hafnium");
    //===========================================================================
    (list[249]).proton = 73;
    (list[249]).neutron = 107;
    (list[249]).mass_of_atom = 179.947466;
    (list[249]).mass_of_nucleus = 327951.0062222962;
    (list[249]).RAM_of_element = 180.9479;
    (list[249]).abundency = 0.00012;
    STR_NAME_COPY((list[249]).symbol, "Ta");
    STR_NAME_COPY((list[249]).name, "Tantalum");
    //===========================================================================
    (list[250]).proton = 73;
    (list[250]).neutron = 108;
    (list[250]).mass_of_atom = 180.947996;
    (list[250]).mass_of_nucleus = 329774.85996272956;
    (list[250]).RAM_of_element = 180.9479;
    (list[250]).abundency = 0.99988;
    STR_NAME_COPY((list[250]).symbol, "Ta");
    STR_NAME_COPY((list[250]).name, "Tantalum");
    //===========================================================================
    (list[251]).proton = 74;
    (list[251]).neutron = 106;
    (list[251]).mass_of_atom = 179.946706;
    (list[251]).mass_of_nucleus = 327948.62082771264;
    (list[251]).RAM_of_element = 183.84;
    (list[251]).abundency = 0.0012;
    STR_NAME_COPY((list[251]).symbol, "W");
    STR_NAME_COPY((list[251]).name, "Tungsten");
    //===========================================================================
    (list[252]).proton = 74;
    (list[252]).neutron = 108;
    (list[252]).mass_of_atom = 181.948206;
    (list[252]).mass_of_nucleus = 331597.13037912763;
    (list[252]).RAM_of_element = 183.84;
    (list[252]).abundency = 0.265;
    STR_NAME_COPY((list[252]).symbol, "W");
    STR_NAME_COPY((list[252]).name, "Tungsten");
    //===========================================================================
    (list[253]).proton = 74;
    (list[253]).neutron = 109;
    (list[253]).mass_of_atom = 182.950224;
    (list[253]).mass_of_nucleus = 333423.69657632464;
    (list[253]).RAM_of_element = 183.84;
    (list[253]).abundency = 0.1431;
    STR_NAME_COPY((list[253]).symbol, "W");
    STR_NAME_COPY((list[253]).name, "Tungsten");
    //===========================================================================
    (list[254]).proton = 74;
    (list[254]).neutron = 110;
    (list[254]).mass_of_atom = 183.950933;
    (list[254]).mass_of_nucleus = 335247.8766136401;
    (list[254]).RAM_of_element = 183.84;
    (list[254]).abundency = 0.3064;
    STR_NAME_COPY((list[254]).symbol, "W");
    STR_NAME_COPY((list[254]).name, "Tungsten");
    //===========================================================================
    (list[255]).proton = 74;
    (list[255]).neutron = 112;
    (list[255]).mass_of_atom = 185.954362;
    (list[255]).mass_of_nucleus = 338899.9025152548;
    (list[255]).RAM_of_element = 183.84;
    (list[255]).abundency = 0.2843;
    STR_NAME_COPY((list[255]).symbol, "W");
    STR_NAME_COPY((list[255]).name, "Tungsten");
    //===========================================================================
    (list[256]).proton = 75;
    (list[256]).neutron = 110;
    (list[256]).mass_of_atom = 184.952956;
    (list[256]).mass_of_nucleus = 337073.45192527515;
    (list[256]).RAM_of_element = 186.207;
    (list[256]).abundency = 0.374;
    STR_NAME_COPY((list[256]).symbol, "Re");
    STR_NAME_COPY((list[256]).name, "Rhenium");
    //===========================================================================
    (list[257]).proton = 75;
    (list[257]).neutron = 112;
    (list[257]).mass_of_atom = 186.955751;
    (list[257]).mass_of_nucleus = 340724.3221161451;
    (list[257]).RAM_of_element = 186.207;
    (list[257]).abundency = 0.626;
    STR_NAME_COPY((list[257]).symbol, "Re");
    STR_NAME_COPY((list[257]).name, "Rhenium");
    //===========================================================================
    (list[258]).proton = 76;
    (list[258]).neutron = 108;
    (list[258]).mass_of_atom = 183.952491;
    (list[258]).mass_of_nucleus = 335248.7166725365;
    (list[258]).RAM_of_element = 190.23;
    (list[258]).abundency = 0.0002;
    STR_NAME_COPY((list[258]).symbol, "Os");
    STR_NAME_COPY((list[258]).name, "Osmium");
    //===========================================================================
    (list[259]).proton = 76;
    (list[259]).neutron = 110;
    (list[259]).mass_of_atom = 185.953838;
    (list[259]).mass_of_nucleus = 338896.94732214714;
    (list[259]).RAM_of_element = 190.23;
    (list[259]).abundency = 0.0159;
    STR_NAME_COPY((list[259]).symbol, "Os");
    STR_NAME_COPY((list[259]).name, "Osmium");
    //===========================================================================
    (list[260]).proton = 76;
    (list[260]).neutron = 111;
    (list[260]).mass_of_atom = 186.955748;
    (list[260]).mass_of_nucleus = 340723.3166474823;
    (list[260]).RAM_of_element = 190.23;
    (list[260]).abundency = 0.0196;
    STR_NAME_COPY((list[260]).symbol, "Os");
    STR_NAME_COPY((list[260]).name, "Osmium");
    //===========================================================================
    (list[261]).proton = 76;
    (list[261]).neutron = 112;
    (list[261]).mass_of_atom = 187.955836;
    (list[261]).mass_of_nucleus = 342546.36467159196;
    (list[261]).RAM_of_element = 190.23;
    (list[261]).abundency = 0.1324;
    STR_NAME_COPY((list[261]).symbol, "Os");
    STR_NAME_COPY((list[261]).name, "Osmium");
    //===========================================================================
    (list[262]).proton = 76;
    (list[262]).neutron = 113;
    (list[262]).mass_of_atom = 188.958145;
    (list[262]).mass_of_nucleus = 344373.4613290834;
    (list[262]).RAM_of_element = 190.23;
    (list[262]).abundency = 0.16149999999999998;
    STR_NAME_COPY((list[262]).symbol, "Os");
    STR_NAME_COPY((list[262]).name, "Osmium");
    //===========================================================================
    (list[263]).proton = 76;
    (list[263]).neutron = 114;
    (list[263]).mass_of_atom = 189.958445;
    (list[263]).mass_of_nucleus = 346196.89580536645;
    (list[263]).RAM_of_element = 190.23;
    (list[263]).abundency = 0.2626;
    STR_NAME_COPY((list[263]).symbol, "Os");
    STR_NAME_COPY((list[263]).name, "Osmium");
    //===========================================================================
    (list[264]).proton = 76;
    (list[264]).neutron = 116;
    (list[264]).mass_of_atom = 191.961479;
    (list[264]).mass_of_nucleus = 349848.20166637516;
    (list[264]).RAM_of_element = 190.23;
    (list[264]).abundency = 0.4078;
    STR_NAME_COPY((list[264]).symbol, "Os");
    STR_NAME_COPY((list[264]).name, "Osmium");
    //===========================================================================
    (list[265]).proton = 77;
    (list[265]).neutron = 114;
    (list[265]).mass_of_atom = 190.960591;
    (list[265]).mass_of_nucleus = 348022.6953321775;
    (list[265]).RAM_of_element = 196.9665;
    (list[265]).abundency = 0.373;
    STR_NAME_COPY((list[265]).symbol, "Ir");
    STR_NAME_COPY((list[265]).name, "Iridium");
    //===========================================================================
    (list[266]).proton = 77;
    (list[266]).neutron = 116;
    (list[266]).mass_of_atom = 192.962924;
    (list[266]).mass_of_nucleus = 351672.72334897163;
    (list[266]).RAM_of_element = 196.9665;
    (list[266]).abundency = 0.627;
    STR_NAME_COPY((list[266]).symbol, "Ir");
    STR_NAME_COPY((list[266]).name, "Iridium");
    //===========================================================================
    (list[267]).proton = 78;
    (list[267]).neutron = 112;
    (list[267]).mass_of_atom = 189.959930;
    (list[267]).mass_of_nucleus = 346197.60279346735;
    (list[267]).RAM_of_element = 192.217;
    (list[267]).abundency = 0.00014000000000000001;
    STR_NAME_COPY((list[267]).symbol, "Pt");
    STR_NAME_COPY((list[267]).name, "Platinum");
    //===========================================================================
    (list[268]).proton = 78;
    (list[268]).neutron = 114;
    (list[268]).mass_of_atom = 191.961035;
    (list[268]).mass_of_nucleus = 349845.39230427635;
    (list[268]).RAM_of_element = 192.217;
    (list[268]).abundency = 0.00782;
    STR_NAME_COPY((list[268]).symbol, "Pt");
    STR_NAME_COPY((list[268]).name, "Platinum");
    //===========================================================================
    (list[269]).proton = 78;
    (list[269]).neutron = 116;
    (list[269]).mass_of_atom = 193.962664;
    (list[269]).mass_of_nucleus = 353494.137008193;
    (list[269]).RAM_of_element = 192.217;
    (list[269]).abundency = 0.32966999999999996;
    STR_NAME_COPY((list[269]).symbol, "Pt");
    STR_NAME_COPY((list[269]).name, "Platinum");
    //===========================================================================
    (list[270]).proton = 78;
    (list[270]).neutron = 117;
    (list[270]).mass_of_atom = 194.964774;
    (list[270]).mass_of_nucleus = 355320.8709110501;
    (list[270]).RAM_of_element = 192.217;
    (list[270]).abundency = 0.33832;
    STR_NAME_COPY((list[270]).symbol, "Pt");
    STR_NAME_COPY((list[270]).name, "Platinum");
    //===========================================================================
    (list[271]).proton = 78;
    (list[271]).neutron = 118;
    (list[271]).mass_of_atom = 195.964935;
    (list[271]).mass_of_nucleus = 357144.05200595537;
    (list[271]).RAM_of_element = 192.217;
    (list[271]).abundency = 0.25242000000000003;
    STR_NAME_COPY((list[271]).symbol, "Pt");
    STR_NAME_COPY((list[271]).name, "Platinum");
    //===========================================================================
    (list[272]).proton = 78;
    (list[272]).neutron = 120;
    (list[272]).mass_of_atom = 197.967876;
    (list[272]).mass_of_nucleus = 360795.18833841634;
    (list[272]).RAM_of_element = 192.217;
    (list[272]).abundency = 0.07163;
    STR_NAME_COPY((list[272]).symbol, "Pt");
    STR_NAME_COPY((list[272]).name, "Platinum");
    //===========================================================================
    (list[273]).proton = 79;
    (list[273]).neutron = 118;
    (list[273]).mass_of_atom = 196.966552;
    (list[273]).mass_of_nucleus = 358968.88722522074;
    (list[273]).RAM_of_element = 195.078;
    (list[273]).abundency = 1.0;
    STR_NAME_COPY((list[273]).symbol, "Au");
    STR_NAME_COPY((list[273]).name, "Gold");
    //===========================================================================
    (list[274]).proton = 80;
    (list[274]).neutron = 116;
    (list[274]).mass_of_atom = 195.965815;
    (list[274]).mass_of_nucleus = 357143.6561470521;
    (list[274]).RAM_of_element = 200.59;
    (list[274]).abundency = 0.0015;
    STR_NAME_COPY((list[274]).symbol, "Hg");
    STR_NAME_COPY((list[274]).name, "Mercury");
    //===========================================================================
    (list[275]).proton = 80;
    (list[275]).neutron = 118;
    (list[275]).mass_of_atom = 197.966752;
    (list[275]).mass_of_nucleus = 360791.13941274275;
    (list[275]).RAM_of_element = 200.59;
    (list[275]).abundency = 0.09970000000000001;
    STR_NAME_COPY((list[275]).symbol, "Hg");
    STR_NAME_COPY((list[275]).name, "Mercury");
    //===========================================================================
    (list[276]).proton = 80;
    (list[276]).neutron = 119;
    (list[276]).mass_of_atom = 198.968262;
    (list[276]).mass_of_nucleus = 362616.7795830338;
    (list[276]).RAM_of_element = 200.59;
    (list[276]).abundency = 0.16870000000000002;
    STR_NAME_COPY((list[276]).symbol, "Hg");
    STR_NAME_COPY((list[276]).name, "Mercury");
    //===========================================================================
    (list[277]).proton = 80;
    (list[277]).neutron = 120;
    (list[277]).mass_of_atom = 199.968309;
    (list[277]).mass_of_nucleus = 364439.7528687515;
    (list[277]).RAM_of_element = 200.59;
    (list[277]).abundency = 0.231;
    STR_NAME_COPY((list[277]).symbol, "Hg");
    STR_NAME_COPY((list[277]).name, "Mercury");
    //===========================================================================
    (list[278]).proton = 80;
    (list[278]).neutron = 121;
    (list[278]).mass_of_atom = 200.970285;
    (list[278]).mass_of_nucleus = 366266.24250466883;
    (list[278]).RAM_of_element = 200.59;
    (list[278]).abundency = 0.1318;
    STR_NAME_COPY((list[278]).symbol, "Hg");
    STR_NAME_COPY((list[278]).name, "Mercury");
    //===========================================================================
    (list[279]).proton = 80;
    (list[279]).neutron = 122;
    (list[279]).mass_of_atom = 201.970626;
    (list[279]).mass_of_nucleus = 368089.7517193439;
    (list[279]).RAM_of_element = 200.59;
    (list[279]).abundency = 0.2986;
    STR_NAME_COPY((list[279]).symbol, "Hg");
    STR_NAME_COPY((list[279]).name, "Mercury");
    //===========================================================================
    (list[280]).proton = 80;
    (list[280]).neutron = 124;
    (list[280]).mass_of_atom = 203.973476;
    (list[280]).mass_of_nucleus = 371740.7221690324;
    (list[280]).RAM_of_element = 200.59;
    (list[280]).abundency = 0.0687;
    STR_NAME_COPY((list[280]).symbol, "Hg");
    STR_NAME_COPY((list[280]).name, "Mercury");
    //===========================================================================
    (list[281]).proton = 81;
    (list[281]).neutron = 122;
    (list[281]).mass_of_atom = 202.972329;
    (list[281]).mass_of_nucleus = 369914.7437069437;
    (list[281]).RAM_of_element = 204.3833;
    (list[281]).abundency = 0.29524;
    STR_NAME_COPY((list[281]).symbol, "Tl");
    STR_NAME_COPY((list[281]).name, "Thallium");
    //===========================================================================
    (list[282]).proton = 81;
    (list[282]).neutron = 124;
    (list[282]).mass_of_atom = 204.974412;
    (list[282]).mass_of_nucleus = 373564.31600183534;
    (list[282]).RAM_of_element = 204.3833;
    (list[282]).abundency = 0.7047599999999999;
    STR_NAME_COPY((list[282]).symbol, "Tl");
    STR_NAME_COPY((list[282]).name, "Thallium");
    //===========================================================================
    (list[283]).proton = 82;
    (list[283]).neutron = 122;
    (list[283]).mass_of_atom = 203.973029;
    (list[283]).mass_of_nucleus = 371737.9073382707;
    (list[283]).RAM_of_element = 207.2;
    (list[283]).abundency = 0.013999999999999999;
    STR_NAME_COPY((list[283]).symbol, "Pb");
    STR_NAME_COPY((list[283]).name, "Lead");
    //===========================================================================
    (list[284]).proton = 82;
    (list[284]).neutron = 124;
    (list[284]).mass_of_atom = 205.974449;
    (list[284]).mass_of_nucleus = 375386.2710586769;
    (list[284]).RAM_of_element = 207.2;
    (list[284]).abundency = 0.24100000000000002;
    STR_NAME_COPY((list[284]).symbol, "Pb");
    STR_NAME_COPY((list[284]).name, "Lead");
    //===========================================================================
    (list[285]).proton = 82;
    (list[285]).neutron = 125;
    (list[285]).mass_of_atom = 206.975881;
    (list[285]).mass_of_nucleus = 377211.7690437344;
    (list[285]).RAM_of_element = 207.2;
    (list[285]).abundency = 0.221;
    STR_NAME_COPY((list[285]).symbol, "Pb");
    STR_NAME_COPY((list[285]).name, "Lead");
    //===========================================================================
    (list[286]).proton = 82;
    (list[286]).neutron = 126;
    (list[286]).mass_of_atom = 207.976636;
    (list[286]).mass_of_nucleus = 379036.03293388;
    (list[286]).RAM_of_element = 207.2;
    (list[286]).abundency = 0.524;
    STR_NAME_COPY((list[286]).symbol, "Pb");
    STR_NAME_COPY((list[286]).name, "Lead");
    //===========================================================================
    (list[287]).proton = 83;
    (list[287]).neutron = 126;
    (list[287]).mass_of_atom = 208.980383;
    (list[287]).mass_of_nucleus = 380864.7509037546;
    (list[287]).RAM_of_element = 208.9804;
    (list[287]).abundency = 1.0;
    STR_NAME_COPY((list[287]).symbol, "Bi");
    STR_NAME_COPY((list[287]).name, "Bismuth");
    //===========================================================================
    (list[288]).proton = 84;
    (list[288]).neutron = 125;
    (list[288]).mass_of_atom = 208.982416;
    (list[288]).mass_of_nucleus = 380867.4568342658;
    (list[288]).RAM_of_element = 209;
    (list[288]).abundency = 0.0;
    STR_NAME_COPY((list[288]).symbol, "Po");
    STR_NAME_COPY((list[288]).name, "Polonium");
    //===========================================================================
    (list[289]).proton = 85;
    (list[289]).neutron = 125;
    (list[289]).mass_of_atom = 209.987131;
    (list[289]).mass_of_nucleus = 382697.9393593469;
    (list[289]).RAM_of_element = 210;
    (list[289]).abundency = 0.0;
    STR_NAME_COPY((list[289]).symbol, "At");
    STR_NAME_COPY((list[289]).name, "Astatine");
    //===========================================================================
    (list[290]).proton = 86;
    (list[290]).neutron = 136;
    (list[290]).mass_of_atom = 222.017570;
    (list[290]).mass_of_nucleus = 404627.0775553077;
    (list[290]).RAM_of_element = 222;
    (list[290]).abundency = 0.0;
    STR_NAME_COPY((list[290]).symbol, "Rn");
    STR_NAME_COPY((list[290]).name, "Radon");
    //===========================================================================
    (list[291]).proton = 87;
    (list[291]).neutron = 136;
    (list[291]).mass_of_atom = 223.019731;
    (list[291]).mass_of_nucleus = 406452.9044254329;
    (list[291]).RAM_of_element = 223;
    (list[291]).abundency = 0.0;
    STR_NAME_COPY((list[291]).symbol, "Fr");
    STR_NAME_COPY((list[291]).name, "Francium");
    //===========================================================================
    (list[292]).proton = 88;
    (list[292]).neutron = 138;
    (list[292]).mass_of_atom = 226.025403;
    (list[292]).mass_of_nucleus = 411930.90667395684;
    (list[292]).RAM_of_element = 226;
    (list[292]).abundency = 0.0;
    STR_NAME_COPY((list[292]).symbol, "Ra");
    STR_NAME_COPY((list[292]).name, "Radium");
    //===========================================================================
    (list[293]).proton = 89;
    (list[293]).neutron = 138;
    (list[293]).mass_of_atom = 227.027747;
    (list[293]).mass_of_nucleus = 413757.06713251467;
    (list[293]).RAM_of_element = 227;
    (list[293]).abundency = 0.0;
    STR_NAME_COPY((list[293]).symbol, "Ac");
    STR_NAME_COPY((list[293]).name, "Actinium");
    //===========================================================================
    (list[294]).proton = 90;
    (list[294]).neutron = 142;
    (list[294]).mass_of_atom = 232.038050;
    (list[294]).mass_of_nucleus = 422889.28639356047;
    (list[294]).RAM_of_element = 232.0381;
    (list[294]).abundency = 1.0;
    STR_NAME_COPY((list[294]).symbol, "Th");
    STR_NAME_COPY((list[294]).name, "Thorium");
    //===========================================================================
    (list[295]).proton = 91;
    (list[295]).neutron = 140;
    (list[295]).mass_of_atom = 231.035879;
    (list[295]).mass_of_nucleus = 421061.44129455916;
    (list[295]).RAM_of_element = 231.0359;
    (list[295]).abundency = 1.0;
    STR_NAME_COPY((list[295]).symbol, "Pa");
    STR_NAME_COPY((list[295]).name, "Protactinium");
    //===========================================================================
    (list[296]).proton = 92;
    (list[296]).neutron = 142;
    (list[296]).mass_of_atom = 234.040946;
    (list[296]).mass_of_nucleus = 426538.340696079;
    (list[296]).RAM_of_element = 238.0289;
    (list[296]).abundency = 5.4999999999999995e-05;
    STR_NAME_COPY((list[296]).symbol, "U");
    STR_NAME_COPY((list[296]).name, "Uranium");
    //===========================================================================
    (list[297]).proton = 92;
    (list[297]).neutron = 143;
    (list[297]).mass_of_atom = 235.043923;
    (list[297]).mass_of_nucleus = 428366.65504249407;
    (list[297]).RAM_of_element = 238.0289;
    (list[297]).abundency = 0.0072;
    STR_NAME_COPY((list[297]).symbol, "U");
    STR_NAME_COPY((list[297]).name, "Uranium");
    //===========================================================================
    (list[298]).proton = 92;
    (list[298]).neutron = 146;
    (list[298]).mass_of_atom = 238.050783;
    (list[298]).mass_of_nucleus = 433847.8228814986;
    (list[298]).RAM_of_element = 238.0289;
    (list[298]).abundency = 0.992745;
    STR_NAME_COPY((list[298]).symbol, "U");
    STR_NAME_COPY((list[298]).name, "Uranium");
    //===========================================================================
    (list[299]).proton = 93;
    (list[299]).neutron = 144;
    (list[299]).mass_of_atom = 237.048167;
    (list[299]).mass_of_nucleus = 432019.16659751086;
    (list[299]).RAM_of_element = 237;
    (list[299]).abundency = 0.0;
    STR_NAME_COPY((list[299]).symbol, "Np");
    STR_NAME_COPY((list[299]).name, "Neptunium");
    //===========================================================================
    (list[300]).proton = 94;
    (list[300]).neutron = 150;
    (list[300]).mass_of_atom = 244.064198;
    (list[300]).mass_of_nucleus = 444807.6025787868;
    (list[300]).RAM_of_element = 244;
    (list[300]).abundency = 0.0;
    STR_NAME_COPY((list[300]).symbol, "Pu");
    STR_NAME_COPY((list[300]).name, "Plutonium");
    //===========================================================================
    (list[301]).proton = 95;
    (list[301]).neutron = 148;
    (list[301]).mass_of_atom = 243.061373;
    (list[301]).mass_of_nucleus = 442978.5653112885;
    (list[301]).RAM_of_element = 243;
    (list[301]).abundency = 0.0;
    STR_NAME_COPY((list[301]).symbol, "Am");
    STR_NAME_COPY((list[301]).name, "Americium");
    //===========================================================================
    (list[302]).proton = 96;
    (list[302]).neutron = 151;
    (list[302]).mass_of_atom = 247.070347;
    (list[302]).mass_of_nucleus = 450285.4743447007;
    (list[302]).RAM_of_element = 247;
    (list[302]).abundency = 0.0;
    STR_NAME_COPY((list[302]).symbol, "Cm");
    STR_NAME_COPY((list[302]).name, "Curium");
    //===========================================================================
    (list[303]).proton = 97;
    (list[303]).neutron = 150;
    (list[303]).mass_of_atom = 247.070299;
    (list[303]).mass_of_nucleus = 450284.3868460954;
    (list[303]).RAM_of_element = 247;
    (list[303]).abundency = 0.0;
    STR_NAME_COPY((list[303]).symbol, "Bk");
    STR_NAME_COPY((list[303]).name, "Berkelium");
    //===========================================================================
    (list[304]).proton = 98;
    (list[304]).neutron = 153;
    (list[304]).mass_of_atom = 251.079580;
    (list[304]).mass_of_nucleus = 457591.8555060038;
    (list[304]).RAM_of_element = 251;
    (list[304]).abundency = 0.0;
    STR_NAME_COPY((list[304]).symbol, "Cf");
    STR_NAME_COPY((list[304]).name, "Californium");
    //===========================================================================
    (list[305]).proton = 99;
    (list[305]).neutron = 153;
    (list[305]).mass_of_atom = 252.082972;
    (list[305]).mass_of_nucleus = 459419.9263507769;
    (list[305]).RAM_of_element = 252;
    (list[305]).abundency = 0.0;
    STR_NAME_COPY((list[305]).symbol, "Es");
    STR_NAME_COPY((list[305]).name, "Einsteinium");
    //===========================================================================
    (list[306]).proton = 100;
    (list[306]).neutron = 157;
    (list[306]).mass_of_atom = 257.095099;
    (list[306]).mass_of_nucleus = 468555.4705588234;
    (list[306]).RAM_of_element = 257;
    (list[306]).abundency = 0.0;
    STR_NAME_COPY((list[306]).symbol, "Fm");
    STR_NAME_COPY((list[306]).name, "Fermium");
    //===========================================================================
    (list[307]).proton = 101;
    (list[307]).neutron = 157;
    (list[307]).mass_of_atom = 258.098425;
    (list[307]).mass_of_nucleus = 470383.4210930143;
    (list[307]).RAM_of_element = 258;
    (list[307]).abundency = 0.0;
    STR_NAME_COPY((list[307]).symbol, "Md");
    STR_NAME_COPY((list[307]).name, "Mendelevium");
    //===========================================================================
    (list[308]).proton = 102;
    (list[308]).neutron = 157;
    (list[308]).mass_of_atom = 259.101024;
    (list[308]).mass_of_nucleus = 472210.04638791265;
    (list[308]).RAM_of_element = 259;
    (list[308]).abundency = 0.0;
    STR_NAME_COPY((list[308]).symbol, "No");
    STR_NAME_COPY((list[308]).name, "Nobelium");
    //===========================================================================
    (list[309]).proton = 103;
    (list[309]).neutron = 159;
    (list[309]).mass_of_atom = 262.109692;
    (list[309]).mass_of_nucleus = 477693.51000771613;
    (list[309]).RAM_of_element = 262;
    (list[309]).abundency = 0.0;
    STR_NAME_COPY((list[309]).symbol, "Lr");
    STR_NAME_COPY((list[309]).name, "Lawrencium");
    //===========================================================================
    (list[310]).proton = 104;
    (list[310]).neutron = 159;
    (list[310]).mass_of_atom = 263.118313;
    (list[310]).mass_of_nucleus = 479531.1127318019;
    (list[310]).RAM_of_element = 261;
    (list[310]).abundency = 0.0;
    STR_NAME_COPY((list[310]).symbol, "Rf");
    STR_NAME_COPY((list[310]).name, "Rutherfordium");
    //===========================================================================
    (list[311]).proton = 105;
    (list[311]).neutron = 157;
    (list[311]).mass_of_atom = 262.011437;
    (list[311]).mass_of_nucleus = 477512.4021855956;
    (list[311]).RAM_of_element = 262;
    (list[311]).abundency = 0.0;
    STR_NAME_COPY((list[311]).symbol, "Db");
    STR_NAME_COPY((list[311]).name, "Dubnium");
    //===========================================================================
    (list[312]).proton = 106;
    (list[312]).neutron = 160;
    (list[312]).mass_of_atom = 266.012238;
    (list[312]).mass_of_nucleus = 484804.4127585712;
    (list[312]).RAM_of_element = 266;
    (list[312]).abundency = 0.0;
    STR_NAME_COPY((list[312]).symbol, "Sg");
    STR_NAME_COPY((list[312]).name, "Seaborgium");
    //===========================================================================
    (list[313]).proton = 107;
    (list[313]).neutron = 157;
    (list[313]).mass_of_atom = 264.012496;
    (list[313]).mass_of_nucleus = 481158.10784357454;
    (list[313]).RAM_of_element = 264;
    (list[313]).abundency = 0.0;
    STR_NAME_COPY((list[313]).symbol, "Bh");
    STR_NAME_COPY((list[313]).name, "Bohrium");
    //===========================================================================
    (list[314]).proton = 108;
    (list[314]).neutron = 161;
    (list[314]).mass_of_atom = 269.001341;
    (list[314]).mass_of_nucleus = 490251.21158228506;
    (list[314]).RAM_of_element = 277;
    (list[314]).abundency = 0.0;
    STR_NAME_COPY((list[314]).symbol, "Hs");
    STR_NAME_COPY((list[314]).name, "Hassium");
    //===========================================================================
    (list[315]).proton = 109;
    (list[315]).neutron = 159;
    (list[315]).mass_of_atom = 268.001388;
    (list[315]).mass_of_nucleus = 488427.4096480027;
    (list[315]).RAM_of_element = 268;
    (list[315]).abundency = 0.0;
    STR_NAME_COPY((list[315]).symbol, "Mt");
    STR_NAME_COPY((list[315]).name, "Meitnerium");
    //===========================================================================
    (list[316]).proton = 110;
    (list[316]).neutron = 162;
    (list[316]).mass_of_atom = 272.001463;
    (list[316]).mass_of_nucleus = 495718.09680457343;
    (list[316]).RAM_of_element = 261.9;
    (list[316]).abundency = 0.0;
    STR_NAME_COPY((list[316]).symbol, "Uun");
    STR_NAME_COPY((list[316]).name, "Darmstadtium");
    //===========================================================================
    (list[317]).proton = 111;
    (list[317]).neutron = 161;
    (list[317]).mass_of_atom = 272.001535;
    (list[317]).mass_of_nucleus = 495717.2280524813;
    (list[317]).RAM_of_element = 271.8;
    (list[317]).abundency = 0.0;
    STR_NAME_COPY((list[317]).symbol, "Uuu");
    STR_NAME_COPY((list[317]).name, "Roentgenium");
    //===========================================================================
    (list[318]).proton = 112;
    (list[318]).neutron = 165;
    (list[318]).mass_of_atom = 277;
    (list[318]).mass_of_nucleus = 504827.86797;
    (list[318]).RAM_of_element = 285;
    (list[318]).abundency = 0.0;
    STR_NAME_COPY((list[318]).symbol, "Uub");
    STR_NAME_COPY((list[318]).name, "Copernicium");
    //===========================================================================
    (list[319]).proton = 114;
    (list[319]).neutron = 175;
    (list[319]).mass_of_atom = 289;
    (list[319]).mass_of_nucleus = 526700.51929;
    (list[319]).RAM_of_element = 289;
    (list[319]).abundency = 0.0;
    STR_NAME_COPY((list[319]).symbol, "Uuq");
    STR_NAME_COPY((list[319]).name, "Flerovium");
    //===========================================================================
    (list[320]).proton = 116;
    (list[320]).neutron = 173;
    (list[320]).mass_of_atom = 289;
    (list[320]).mass_of_nucleus = 526698.51929;
    (list[320]).RAM_of_element = 293;
    (list[320]).abundency = 0.0;
    STR_NAME_COPY((list[320]).symbol, "Uuh");
    STR_NAME_COPY((list[320]).name, "Livermorium");
    //===========================================================================
    (list[321]).proton = 118;
    (list[321]).neutron = 175;
    (list[321]).mass_of_atom = 293;
    (list[321]).mass_of_nucleus = 533988.06973;
    (list[321]).RAM_of_element = 294;
    (list[321]).abundency = 0.0;
    STR_NAME_COPY((list[321]).symbol, "Uuo");
    STR_NAME_COPY((list[321]).name, "Ununoctium");
    //===========================================================================

    return;
}
double AU_UNIT_MASS(int amount_proton, int amount_neutron)
{
    double mass;
    struct mass_list list[_RAM_LIST_LENGTH_];
    int i;
    mass = -1.0;
    amount_proton = (amount_proton >= 0 ? amount_proton : -amount_proton);
    amount_neutron = (amount_neutron >= 0 ? amount_neutron : -amount_neutron);
    if (list_if_inilialized == 0)
    {
        value_mass_list(list);
        list_if_inilialized = 1;
    }
    if (amount_proton == 0)
        mass = (double)amount_neutron * (list[0]).mass_of_nucleus;
    else
        for (i = 1; i < _RAM_LIST_LENGTH_; i++)
        if ((list[i]).proton == amount_proton && (list[i]).neutron == amount_neutron)
            mass = (list[i]).mass_of_nucleus;
    if (mass < 0.0)
    {
        mass = (double)(amount_proton + amount_neutron) * 0.25 * (list[5]).mass_of_nucleus;
        mass = mass - (double)amount_proton;
    }
    return mass; 
}
static void STR_NAME_COPY(char*to, char*from)
{
    int i;
    for (i = 0; from[i] != '\0'; i++)
        to[i] = from[i];
    to[i] = '\0';
    return;
}
//=====================================================================================
//-----------------------------------wasted code---------------------------------------
//=====================================================================================
//--------------------------------------No.1-------------------------------------------
/*
//==================================Hydrogen=================================
    (list[1]).proton = 1;
    (list[1]).neutron = 0;
    (list[1]).mass = 1.007825;//proton 1.00727642
    //=================================Deuterium=================================
    (list[2]).proton = 1;
    (list[2]).neutron = 1;
    (list[2]).mass = 2.014102;
    //==================================Tritium==================================
    (list[3]).proton = 1;
    (list[3]).neutron = 2;
    (list[3]).mass = 3.016049;
    //====================================3He====================================
    (list[4]).proton = 2;
    (list[4]).neutron = 1;
    (list[4]).mass = 3.016029;
    //====================================4He====================================
    (list[5]).proton = 2;
    (list[5]).neutron = 2;
    (list[5]).mass = 4.002603;
    //====================================6Li====================================
    (list[6]).proton = 3;
    (list[6]).neutron = 3;
    (list[6]).mass = 6.015122;
    //====================================7Li====================================
    (list[7]).proton = 3;
    (list[7]).neutron = 4;
    (list[7]).mass = 7.016004;
    //====================================9Be====================================
    (list[8]).proton = 4;
    (list[8]).neutron = 5;
    (list[8]).mass = 9.012182;
    //====================================10B====================================
    (list[9]).proton = 5;
    (list[9]).neutron = 5;
    (list[9]).mass = 10.012937;
    //====================================11B====================================
    (list[10]).proton = 5;
    (list[10]).neutron = 6;
    (list[10]).mass = 11.009305;
    //====================================12C====================================
    (list[11]).proton = 6;
    (list[11]).neutron = 6;
    (list[11]).mass = 12.000000;
    //====================================13C====================================
    (list[12]).proton = 6;
    (list[12]).neutron = 7;
    (list[12]).mass = 13.003355;
    //====================================14C====================================
    (list[13]).proton = 6;
    (list[13]).neutron = 8;
    (list[13]).mass = 14.003242;
    //====================================14N====================================
    (list[14]).proton = 7;
    (list[14]).neutron = 7;
    (list[14]).mass = 14.003074;
    //====================================15N====================================
    (list[15]).proton = 7;
    (list[15]).neutron = 8;
    (list[15]).mass = 15.000109;
*/