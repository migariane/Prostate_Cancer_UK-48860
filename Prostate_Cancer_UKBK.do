/* Metabolic syndrome paradox in prostate cancer: A UK Biobank prospective cohort study
=================================================================================================
Shing Fung Lee 1 2, Maja Nikšić 3,4 , Miguel Angel Luque-Fernandez 3, 5 (Senior and corresponding author)
Authors' affiliations
1 Department of Radiation Oncology, National University Cancer Institute, National University Hospital, Singapore
2 Department of Clinical Oncology, Tuen Mun Hospital, New Territories West Cluster, Hospital Authority, Hong Kong
3 Department of Non-Communicable Disease Epidemiology, London School of Hygiene and Tropical Medicine, London, United Kingdom
4 Centre for Health Services Studies, University of Kent, Canterbury, United Kingdom
5 Department of Statistics and Operations Research, University of Granada, Granada, Spain 

Copyright (c) 2023 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Bug reports: miguel-angel.luque at lshtm.ac.uk	
*/

/////////////////////////////////////////////////////////////////////////////
// DATA MANAGEMENT AND ANALYSIS for Prostate cancer UKBK project # 48860
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// DATA MANAGEMENT 
/////////////////////////////////////////////////////////////////////////////

set processor 16

//////////////////////////////////////////////
// WE MUST DELETE UKBK WITHDRAWS IF ANY
//////////////////////////////////////////////

count
local id 1010995 1016806 1038349 1055033 1157020 1221515 1273024 1275336 1313156 1515177 1664804 1844971 2223588 2378503 2881075 2939084 3119798 3120595 3334692 3335524 3347919 3419048 3426357 3489310 3561476 3737718 3904295 4009646 4107243 4110030 4303119 4501877 4527447 4742847 4759280 4801062 5177119 5314963 5416920 5529277 5733201 5810242 5827597 5907535 5973661
foreach i of local id {
	drop if eid==`i'
				 }
count				 
 // In our data file, we have to exclude participants who withdrew from the study

* (MN) Total N=109,441 (difference is 22 patients)

// To include men only
gen sex =  n_31_0_0    // Sex
keep if sex==1   //keep male participants
count // N= 50,491


// Cancer SITES (Probably the most relevant and important variable)
//First cancer
tab s_40006_0_0
icd10 clean s_40006_0_0, generate(cancersite) nodots  // cancer site icd10
icd10 generate csite = cancersite, description  addcode(begin) nodots // Cancer site labels
icd10 generate catcode = cancersite, category
tab catcode // C61 = prostate cancer

gen cancer_date = ts_40005_0_0    // ts_40005_0_0 is DATE of first cancer diagnosis
format cancer_date %td
gen assessment_date = ts_53_0_0   // ts_53_0_0 is DATE of first attending UKBiobank assessment
format assessment_date %td


// We can include patient with first dx of non-melanotic skin cancer, for these patients they are followed until censored or development of prostate cancer (s40006_1_0 = second episode of cancer)
tab s_40006_1_0
icd10 clean s_40006_1_0, generate(cancersite2) nodots  // cancer site icd10
icd10 generate csite2 = cancersite2, description  addcode(begin) nodots // Cancer site labels
icd10 generate catcode2 = cancersite2, category
tab catcode2 if catcode2=="C61" // C61 = prostate cancer

gen cancer_date2 = ts_40005_1_0    // ts_40005_0_0 is DATE of second cancer diagnosis (for those with first cancer as non-melanotic skin cancer)
format cancer_date2 %td

// To exclude prevalence cases (developed cancer before UK Biobank follow-up), so that all cancer cases were diagnosed during UK Biobank follow-up
gen ci = cancer_date - assessment_date if catcode!="C44" //date at 1st cancer diagnosis - date at UK biobank baseline recruitment == INCIDENT CANCER CASES
gen ci2 = cancer_date2 - assessment_date

tab ci if ci>1   // To avoid reverse casuality, we can analyze only the incidence prostate cancer cases. Total incident cancer cases (including prostate and non-prostate cancers = 28,533)
tab ci if ci>=0

drop if ci<0 //keep incident cancer cases & participants who did not develop cancer during followup = 242,524 (all men in UKBiobank with no baseline cancer, except non-melanoma skin cancers [C44])

drop if ci2<0 //drop those with skin melanoma then prostate cancer but both developed before UK biobank.

count    
tab ci, missing  // "missing" = no incident cancer during follow-up

gen id = _n   // generate ID

gen prostate = 0
replace prostate = 1 if catcode=="C61"| catcode2=="C61"
replace prostate = 0 if  n_40012_0_0==0| n_40012_0_0==1 |n_40012_0_0==2   //Exclude benign, Uncertain whether benign or malignant, and carcinoma in situ
replace prostate = 0 if  n_40012_1_0==0| n_40012_1_0==1 |n_40012_1_0==2

gen prostate_date = cancer_date if prostate==1
replace prostate_date = cancer_date2 if catcode=="C44" & catcode2=="C61"& prostate==1 

gen lastfu_prostate = cancer_date if catcode=="C61"
replace lastfu_prostate = cancer_date2 if catcode2=="C61" & lastfu_prostate==.
replace lastfu_prostate =date("2021, 1, 31", "YMD") if lastfu_prostate==.    //  last censoring date on 31 Jan 2021
replace lastfu_prostate =date("2019, 7, 31", "YMD") if lastfu_prostate==.    //  last censoring date on 31 July 2019


gen death = 1 if date_of_death_2!=.
gen DOD = date_of_death_2

gen competing_risk = 0    // 0 = no prostate cancer or other cancer, no death; 1 = prostate cancer; 2 = other cancers or any death
replace competing_risk = 1 if prostate ==1 
replace competing_risk = 2 if death==1 & s_40006_0_0=="" 
replace competing_risk = 2 if death==0 & prostate ==0 & s_40006_0_0!=""

gen competing_risk_date = prostate_date
replace competing_risk_date = DOD if death==1 & s_40006_0_0==""  
replace competing_risk_date = cancer_date if prostate==0 & s_40006_0_0!=""  
format %tdDD/NN/CCYY competing_risk_date

gen competing_risk_date_new = prostate_date  // last censoring date on 31 July 2019
format %tdDD/NN/CCYY competing_risk_date_new
replace competing_risk_date_new = DOD if death==1 & s_40006_0_0=="" & DOD <= date("2019, 7, 31", "YMD")
replace competing_risk_date_new = cancer_date if prostate==0 & s_40006_0_0!=""  & cancer_date <= date("2019, 7, 31", "YMD")
replace competing_risk_date_new = date("2019, 7, 31", "YMD") if competing_risk_date_new ==.





// Quintiles of Townsend deprivation score (qdi)
gen townsend = n_189_0_0
xtile qdi = n_189_0_0, nq(5) // Quintiles Townsend deprivation score 

/**   we will use Townsend deprivation score instead

// Index of Multiple Deprivation IMD (England, Scotland, and Wales)
gen IMD_Eng = n_26410_0_0
gen IMD_Scot = n_26427_0_0
gen IMD_Wales = n_26426_0_0
gen IDT = n_189_0_0

xtile IMD_Eng5 = IMD_Eng, nq(5) // Quintiles of IMD_england
label variable IMD_Eng5 "Socioeconomic deprivation quintiles in IMD_England"
tab IMD_Eng5 ccst, miss col chi
tab IMD_Eng5 ccst, col chi

* Quintiles of IMD_Eng by country
xtile IMD_Eng5 = IMD_Eng, nq(5)
xtile IMD_Wales5 = IMD_Wales, nq(5)
xtile IMD_Scot5 = IMD_Scot, nq(5)

tab IMD_Eng5, miss
tab IMD_Wales5, miss
tab IMD_Scot5, miss

* Labeling IMD 
lab def imd 1"Most affluent" 2"Affluent" 3"Average" 4"Deprived" 5"Most deprived"
label val IMD_Eng5 imd
label val IMD_Wales5 imd
label val IMD_Scot5 imd

**/

// Ethnicity (ethnic)
tab n_21000_0_0
gen ethnic = n_21000_0_0
recode ethnic -3=. -1 =.
recode ethnic 1001/1003=1    2001/2004=2   3001/3004=3  4001/4003=4    5=3 6=5
tab ethnic
label variable ethnic "1= white, 2 = mixed/others, 3 = asian (including Chinese), 4 = black
lab def ethnic 1"white" 2"mixed/others" 3"asian including Chinese" 4"black"
label val ethnic ethnic
recode ethnic 5 = 2 
tab ethnic

// Education qualifications (quali)
tab n_6138_0_0
gen quali =  n_6138_0_0
recode quali -3=. -7 =.   // -7: None of the above; -3: Prefer not to answer
recode quali  6=1 2/5=2  //   dichotomous variable comparing degree level education + professional qualifications with all other qualifications
label def quali 1"college/degree level education + professional qualifications" 2"all other qualifications"
label variable quali "1=college/degree level education + professional qualifications, 2=all other qualifications"

// Job 
tab n_6142_0_0
gen job = n_6142_0_0
recode job -3=. -7 =.   // -7: None of the above; -3: Prefer not to answer
recode job 3/7 =3
tab job
label variable job "1= employed or self-employed, 2 = retired, 3 = unemployed/unpaid/students"
label def job 1"employed or self-employed" 2"retired" 3"unemployed/unpaid/students"

//marital status (marital)    // derived from self-reported household occupancy and relatedness data as follows: married/partner was derived from those reporting husband/wife/partner in household with >1 person reported to live in household
tab n_709_0_0   // no. of household members
gen household =  n_709_0_0
recode household -3=. -1 =.  1=0   2/100=1  

tab n_6141_0_0   // wife/husband/partner = 1
gen relation = n_6141_0_0
recode relation -3=.  2/10=0

gen marital = relation if household ==1
tab marital

// Age
gen age =  2021 - n_34_0_0 // n_34_0_0 is year of birth 
gen age_cat = age
recode age_cat 50/65=1 66/70=2 71/75=3 76/80=4 81/90=5

// BMI 
gen bmi =  n_21001_0_0 // BMI in kg/m2
summarize bmi,detail   
gen bmi_cat = bmi
recode bmi_cat 0/18.4999=0 18.5000/24.9999=1 25/29.99999=2 30/100=3
recode bmi_cat 0/1=1

// smoking status (smoking)   // 0 = non-smoker, 1 = ex-smoker, 2 = current smoker, 3 = prefer not to answer
tab n_20116_0_0
gen smoking = n_20116_0_0
recode smoking -3=. 
lab def smoking  0"non-smoker" 1"ex-smoker" 2"current smoker" 


// International Physical Activity Questionnaire (ipaq)  0 = low, 1 = moderate, 2 = high
tab n_22032_0
gen ipaq = n_22032_0
 ** MN: recode ipad -3=. //An error should be ipaq
recode ipaq -3=.

// Alcohol intake frequecy
tab n_1558_0_0
gen alc_freq = n_1558_0_0
recode alc_freq -3=.

// Portion of fruit and vetegable
tab n_1309_0_0
gen fruit_freq = n_1309_0_0  
recode fruit_freq -3=. -1=. -10=0   // -10 (Less than one), -3 (Prefer not to answer), -1 (Do not know)
recode fruit_freq 0/4= 0 5/100=1 

// Processed meat
tab n_1349_0_0
gen processedmeat_freq = n_1349_0_0  
recode processedmeat_freq -3=. -1=.  // -3 (Prefer not to answer), -1 (Do not know)  0	Never, 1=Less than once a week, 2= Once a week, 3=2-4 times a week, 4= 5-6 times a week, 5 =	Once or more daily

// Level of physical activity (0 = low, 1= moderate, and 2= high)
tab n_22032_0_0
gen physical = n_22032_0_0

// Ever had prostate specific antigen (PSA) test
tab n_2365_0_0       
gen ever_PSA = n_2365_0_0 
recode ever_PSA -3=. -1=.     //-1 Do not know, -3 prefer not to answer

// Family history of prostate cancer
tab n_20107_0_0
gen father_prostate = 0
replace father_prostate = 1 if n_20107_0_0 == 13      // father had prostate cancer
replace father_prostate = . if n_20107_0_0 == .  

tab n_20111_0_0
gen sibling_prostate = 0
replace sibling_prostate = 1 if n_20111_0_0 == 13    // sibling had prostate cancer
replace sibling_prostate = . if n_20111_0_0 == .

// IGF-I
summarize n_30770_0_0, detail
gen IGF = n_30770_0_0
xtile IGF_med= IGF
gen IGF_mean =1
replace IGF_mean =2 if IGF> 21.70388
 
// Testosterone
summarize n_30850_0_0, detail
gen testosterone = n_30850_0_0
xtile testosterone_med= testosterone
gen testosterone_mean =1
replace testosterone_mean =2 if testosterone> 11.98283

//ICD-10 diagnoses
/*
Hypertension: I10-I16	
obesity/overweight: E66, Z68.3, Z68.4,  O99.21, R93.9	
Hyperlipidaemia: E78.0,  E78.1, E78.2, E78.3, E78.4, E78.5	
hyperglycaemia/Diabetes: E08 to E13
*/

tab s_41202_0_0  //main ICD10 diagnosis
icd10 clean s_41202_0_0, generate(main_dx) nodots  
icd10 generate main_dx_cat = main_dx, category

tab s_41204_0_0   //secondary ICD10 diagnosis
icd10 clean s_41204_0_0, generate(secondary_dx) nodots 
icd10 generate secondary_dx_cat = secondary_dx, category 


gen icd_HTN = 0  //hypertension on ICD-10
replace icd_HTN = 1 if main_dx_cat == "I10" |main_dx_cat == "I11" |main_dx_cat == "I12" |main_dx_cat == "I13" |main_dx_cat == "I14" |main_dx_cat == "I15" |main_dx_cat == "I16" 
replace icd_HTN = 1 if secondary_dx_cat == "I10" |secondary_dx_cat == "I11" |secondary_dx_cat == "I12" |secondary_dx_cat == "I13" |secondary_dx_cat == "I14" |secondary_dx == "I15" |secondary_dx_cat == "I16"

gen icd_hyperlipid = 0     // hyperlipidemia on ICD-10
replace icd_hyperlipid = 1 if main_dx == "E780" |main_dx == "E781" |main_dx == "E782" |main_dx == "E783" |main_dx == "E784" |main_dx == "E785" 
replace icd_hyperlipid = 1 if secondary_dx == "E780" |secondary_dx == "E781" |secondary_dx == "E782" |secondary_dx == "E783" |secondary_dx == "E784" |secondary_dx == "E785" 

gen icd_DM = 0    // diabetes on ICD-10
replace icd_DM = 1 if main_dx_cat == "E08" |main_dx_cat == "E09" |main_dx_cat == "E10" |main_dx_cat == "E11" |main_dx_cat == "E12"| main_dx_cat == "E13"
replace icd_DM = 1 if secondary_dx_cat == "E08" |secondary_dx_cat == "E09" |secondary_dx_cat == "E10" |secondary_dx_cat == "E11" |secondary_dx_cat == "E12"| secondary_dx_cat == "E13"

gen icd_obesity = 0   // obesity or overweight on ICD-10
replace icd_obesity = 1 if main_dx_cat == "E66" | main_dx == "Z683"|main_dx == "Z684"| main_dx == "O9921"| main_dx == "R939"
replace icd_obesity = 1 if main_dx_cat == "E66" | secondary_dx == "Z683"|secondary_dx  == "Z684"| secondary_dx == "O9921"| secondary_dx == "R939"


// diagnosis of Metabolic Syndrome by "Modified NCEP ATP III guideline"
/*
Hypertension: systolic blood pressure >=130 and/or diastolic blood pressure >=85
obesity/overweight: waist circumference >= 102
Hyperlipidaemia: Triglyceride >= 1.7, High-density lipoprotein < 1.03
hyperglycaemia/Diabetes: HbA1c >=42
*/
gen SBP = n_4080_0_0
gen SBP_cat = SBP
recode SBP_cat 60/119=0 120/129=1 130/139=2 140/300=3

gen DBP = n_4079_0_0
gen DBP_cat = DBP
recode DBP_cat 0/84.99999=0 85/300=1

gen HDL = n_30760_0_0
gen TG = n_30870_0_0
gen hba1c = n_30750_0_0

gen crp = n_30710_0_0   //C-reactive protein
summarize crp, detail

gen hypertension = 0
replace hypertension =1 if icd_HTN==1|SBP>=130|DBP>=85|medications==2

gen hyperlipidemia = 0
replace hyperlipidemia =1 if icd_hyperlipid==1|TG>=1.7|medications==1

gen low_HDL = 0
replace low_HDL =1 if HDL <1.03 

gen DM = 0
replace DM =1 if icd_DM==1|hba1c >=42|medications==3

gen obese = 0
replace obese =1 if icd_obesity==1|wc>=102

// Use of medications (1 = Cholesterol lowering medication, 2 = Blood pressure medication, 3=Insulin, -7 = none of the above, -1 = Do not know. -3= Prefer not to answer)
tab n_6177_0_0
gen medications = n_6177_0_0
recode medications -3=. -1 =. -7=0  //-1 "Do not know", -3 "Prefer not to answer", -7 "None of the above"    1=Cholesterol lowering medication, 2=Blood pressure medication, 3=Insulin
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// DATA ANALYSIS 
////////////////////////////////////////////////////////////////////////////////

// Table 1 Baseline Characteristics

//Number of metabolic syndrome components
gen MetSyn = hypertension + hyperlipidemia + low_HDL + DM + obese
tab MetSyn, missing
tab MetSyn if prostate ==1, missing
gen MetSyn_new = MetSyn

recode MetSyn_new 4/5=4 //4 or above =4
tab MetSyn_new, missing
tab MetSyn_new if prostate ==1, missing
tab MetSyn_new if prostate ==0, missing
tab MetSyn_new prostate,  col chi

//Metabolic syndrome yes versus no
tab MetSyn_bi, missing
tab MetSyn_bi if prostate ==1, missing
tab MetSyn_bi if prostate ==0, missing
tab MetSyn_bi prostate,  col chi


//Individual metabolic syndrome components
tab hypertension, missing
tab hypertension if prostate ==1, missing
tab hypertension if prostate ==0, missing
tab hypertension prostate,  col chi

tab hyperlipidemia, missing
tab hyperlipidemia if prostate ==1, missing
tab hyperlipidemia if prostate ==0, missing
tab hyperlipidemia prostate,  col chi

tab low_HDL, missing
tab low_HDL if prostate ==1, missing
tab low_HDL if prostate ==0, missing
tab low_HDL prostate,  col chi

tab DM, missing
tab DM if prostate ==1, missing
tab DM if prostate ==0, missing
tab DM prostate,  col chi

tab obese
tab obese if prostate ==1, missing
tab obese if prostate ==0, missing
tab obese prostate,  col chi

//Sociodemographic
summarize age, detail //Age at recruitment (years), mean (SD) 
summarize age if prostate ==1, detail
summarize age if prostate ==0, detail
ttest age, by (prostate)
ranksum age, by(prostate)

tab marital
tab marital if prostate ==1, missing
tab marital if prostate ==0, missing
tab marital prostate,  col chi

tab ethnic
tab ethnic if prostate ==1, missing // 4= Black ethnicity
tab ethnic if prostate ==0, missing 
tab ethnic prostate, col chi

tab qdi
tab qdi if prostate ==1, missing // 5= most deprived, 1= most affluent townsend
tab qdi if prostate ==0, missing 
tab qdi prostate, col chi

tab smoking, missing
tab smoking if prostate ==1, missing //  0 "non-smoker" 1 "ex-smoker" 2 "current smoker" 
tab smoking if prostate ==0, missing 
tab smoking prostate, col chi

recode processedmeat_freq 0/1=1 2/3=2 4/5=3    //0	Never, 1=Less than once a week, 2= Once a week, 3=2-4 times a week, 4= 5-6 times a week, 5 =	Once or more daily

tab processedmeat_freq, missing     
tab processedmeat_freq if prostate ==1, missing
tab processedmeat_freq if prostate ==0, missing
tab processedmeat_freq prostate, col chi

tab fruit_freq, missing
tab fruit_freq if prostate ==1, missing
tab fruit_freq if prostate ==0, missing
tab fruit_freq prostate, col chi

tab ipaq prostate, missing col chi  //physical exercise level: low mod high

summarize bmi, detail
summarize bmi if prostate ==1, detail
summarize bmi if prostate ==0, detail
ttest bmi, by (prostate)

tab medications, missing  // Self-reported medications (1 = Cholesterol lowering medication, 2 = Blood pressure medication, 3=Insulin, -7 = none of the above, -1 = Do not know. -3= Prefer not to answer)
tab medications prostate, missing col chi

tab ever_PSA, missing
tab ever_PSA if prostate==1, missing   // Ever had prostate specific antigen (PSA) test
tab ever_PSA if prostate==0, missing
tab ever_PSA prostate, missing col chi

tab father_prostate, missing
tab father_prostate if prostate==1, missing   // Father had CA prostate
tab father_prostate if prostate==0, missing
tab father_prostate prostate,  col chi

tab sibling_prostate, missing
tab sibling_prostate if prostate==1, missing   // Siblings had CA prostate
tab sibling_prostate if prostate==0, missing
tab sibling_prostate prostate, missing col chi

tab crp if crp <1, missing
tab crp if crp >=1, missing
tab crp if crp <1 & prostate==1, missing
tab crp if crp >=1 & prostate==1, missing

gen crp_cat = crp
recode crp_cat 0/1=0 1.000001/100=1
tab crp_cat 

tab crp_cat, missing
tab crp_cat if prostate ==1, missing
tab crp_cat if prostate ==0, missing
tab crp_cat prostate, missing col chi

tab testosterone_mean prostate, missing col chi

tab IGF_mean prostate, missing col chi





// Median follow-up time from the baseline assessment foir the whole cohort
summarize _t, detail
display in smcl as text "Median follow-up (years): " as result %3.0f r(p50)

// Median interval to prostate cancer incidence from date of baseline assessment 
stsum, by(prostate)


// Incidence rate of prostate cancer by variables

stset lastfu_prostate, failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    

stptime 
/*
Cohort     |  Person-time   Failures        Rate   [95% Conf. Interval]
-----------+-----------------------------------------------------------
     Total |    2894638.8       7038   .00243139   .0023752    .0024889

*/

// gen ethnic_cat = ethnic
// recode ethnic_cat 2/3=2 4=3

local predictorvar "age_cat qdi ethnic marital bmi_cat smoking fruit_freq processedmeat_freq ipaq hypertension hyperlipidemia low_HDL DM obese MetSyn_bi MetSyn_new ever_PSA father_prostate sibling_prostate testosterone_mean IGF_mean"
foreach p of local predictorvar {
stptime, by(`p') per(100000)
}

// Piecewise Exponential Poisson Method (Univariable analysis, for Table 2) to obtain the IRR

local predictorvar "age_cat qdi ethnic marital bmi_cat smoking fruit_freq processedmeat_freq ipaq hypertension hyperlipidemia low_HDL DM obese MetSyn_bi MetSyn_new ever_PSA father_prostate sibling_prostate testosterone_mean IGF_mean"

foreach p of local predictorvar {
preserve 

stsplit split_time, at(1(1)15)
generate risktime = _t - _t0
list id _t0 _t prostate split_time risktime if inlist(id, 1,2), sepby(eid) ab(20) noobs

collapse (min) start=_t0 (max) end=_t (count) n=prostate (sum) risktime prostate, by(split_time `p')

generate midt = (start + end)/2
rcsgen midt, df(3) gen(t_rcs) fw(prostate) orthog    
local knots 'r (knots)'

glm prostate t_rcs1-t_rcs3 i.`p', family (poisson) lnoffset (risktime) nolog noheader eform
restore

}

// To create standardized variables for LASSO
local predictorvar "age_cat qdi ethnic_cat quali job marital bmi_cat wc_cat whratio_cat smoking alc_freq_new_2 fruit_freq ipaq hypertension hyperlipidemia low_HDL DM obese MetSyn_new TG_cat hba1c_cat SBP_cat DBP_cat crp_cat shift_work shift_job physical_work work_asbestos work_paints work_pesticides work_per_week_cat work_sat_new ever_PSA father_prostate sibling_prostate children_no firstsex_cat sexpartner_no homosex homosexpartner_no"

foreach p of local predictorvar {
egen z2_`p' = std(`p')
}

// To run lasso and reduce dimensionality use
preserve 

stsplit split_time, at(1(1)15)
generate risktime = _t - _t0
list id _t0 _t prostate split_time risktime if inlist(id, 1,2), sepby(eid) ab(20) noobs

collapse (min) start=_t0 (max) end=_t (count) n=prostate (sum) risktime prostate, by(split_time z2_age_cat z2_qdi z2_ethnic z2_quali z2_job z2_marital z2_bmi_cat z2_wc_cat z2_whratio_cat z2_smoking z2_alc_freq_new_2 z2_fruit_freq z2_ipaq z2_hypertension z2_hyperlipidemia z2_low_HDL z2_DM z2_obese z2_MetSyn_new z2_TG_cat z2_hba1c_cat z2_SBP_cat z2_DBP_cat z2_crp_cat z2_shift_work z2_shift_job z2_physical_work z2_work_asbestos z2_work_paints z2_work_pesticides z2_work_per_week_cat z2_work_sat_new z2_ever_PSA z2_father_prostate z2_sibling_prostate z2_children_no z2_firstsex_cat z2_sexpartner_no z2_homosex z2_homosexpartner_no)  // ....,by(all the coariates for analysis here…)

generate midt = (start + end)/2
rcsgen midt, df(3) gen(t_rcs) fw(prostate) orthog    
local knots 'r (knots)'


// LASSO		
//split data  /*---------- Step 1: split data --------------*/
splitsample, generate(sample) split(0.9 0.1) rseed(12345)
label define lbsample 1 "Training" 2 "Validation"
label value sample lbsample
tabulate sample

//glm prostate age_cat i.qdi i.ethnic quali household whratio i.smoking fruit_freq i.physical hypertension hyperlipidemia i.MetSyn TG SBP_cat i.medications i.shift_job work_per_week i.work_sat ever_PSA father_prostate sibling_prostate children_no firstsex sexpartner_no  homosexpartner_no, family (poisson) offset (risktime) nolog noheader eform vce(robust)   

lasso poisson prostate z2_hypertension z2_hyperlipidemia z2_low_HDL z2_DM z2_obese z2_MetSyn_new z2_TG_cat z2_hba1c_cat z2_SBP_cat z2_DBP_cat z2_crp_cat t_rcs1-t_rcs3 z2_age_cat z2_qdi z2_ethnic z2_quali z2_job z2_marital z2_bmi_cat z2_wc_cat z2_whratio_cat z2_smoking z2_alc_freq_new_2 z2_fruit_freq z2_ipaq  z2_shift_work z2_shift_job z2_physical_work z2_work_asbestos z2_work_paints z2_work_pesticides z2_work_per_week_cat z2_work_sat_new z2_ever_PSA z2_father_prostate z2_sibling_prostate z2_children_no z2_firstsex_cat z2_sexpartner_no z2_homosex z2_homosexpartner_no if sample==1, exposure(risktime) nolog selection(cv) rseed(12345)     /*---------- Step 2: run in training sample ----*/

cvplot, minmax name(lasso_cvplot, replace) // shows the cross-validation error of different lambda* values
estimates store cv
lassoknots, display (nonzero)

lasso poisson prostate z2_hypertension z2_hyperlipidemia z2_low_HDL z2_DM z2_obese z2_MetSyn_new z2_TG_cat z2_hba1c_cat z2_SBP_cat z2_DBP_cat z2_crp_cat t_rcs1-t_rcs3 z2_age_cat z2_qdi z2_ethnic z2_quali z2_job z2_marital z2_bmi_cat z2_wc_cat z2_whratio_cat z2_smoking z2_alc_freq_new_2 z2_fruit_freq z2_ipaq  z2_shift_work z2_shift_job z2_physical_work z2_work_asbestos z2_work_paints z2_work_pesticides z2_work_per_week_cat z2_work_sat_new z2_ever_PSA z2_father_prostate z2_sibling_prostate z2_children_no z2_firstsex_cat z2_sexpartner_no z2_homosex z2_homosexpartner_no if sample==1, exposure(risktime) nolog selection(adaptive) rseed(12345)  

cvplot, minmax name(lasso_adapativeplot, replace) // shows the cross-validation error of different lambda* values
estimates store adaptive
lassoknots, display (nonzero)

lassoinfo cv adaptive
lassocoef cv adaptive, sort (coef, standardized)   // shows which variables are retained

lassogof cv adaptive, over(sample) postselection   /*---------- Step 3: Evaluate prediciton in testing sample with goodness of fit ----*/

restore


/*
. lassocoef cv adaptive, sort (coef, standardized)   // shows which variables are retained
-------------------------------------------
                     |    cv      adaptive 
---------------------+---------------------
               _cons |     x         x     
              t_rcs1 |     x         x     
              t_rcs2 |     x         x     
z2_homosexpartner_no |     x         x     
          z2_homosex |     x         x     
          z2_age_cat |     x         x     
        z2_hba1c_cat |     x         x     
       z2_shift_work |     x         x     
               z2_DM |     x         x     
              t_rcs3 |     x         x     
          z2_smoking |     x         x     
             z2_ipaq |     x         x     
       z2_MetSyn_new |     x         x     
          z2_SBP_cat |     x         x     
           z2_wc_cat |     x         x     
   z2_alc_freq_new_2 |     x         x     
    z2_physical_work |     x         x     
    z2_work_asbestos |     x         x     
     z2_hypertension |     x         x     
  z2_father_prostate |     x         x     
       z2_ethnic_cat |     x         x     
          z2_crp_cat |     x         x     
          z2_DBP_cat |     x         x     
z2_work_per_week_cat |     x         x     
  z2_work_pesticides |     x         x     
   z2_hyperlipidemia |     x    
          z2_marital |     x         x     
          z2_bmi_cat |     x    
       z2_fruit_freq |     x         x     
 z2_sibling_prostate |     x         x     
      z2_children_no |     x    
      z2_whratio_cat |     x    
        z2_shift_job |     x    
     z2_firstsex_cat |     x    
      z2_work_paints |     x    
           z2_TG_cat |     x    
         z2_ever_PSA |     x    
              z2_qdi |     x    
          z2_low_HDL |     x    
    z2_sexpartner_no |     x    
            z2_quali |     x    
     z2_work_sat_new |     x    
-------------------------------------------
Legend:
  b - base level
  e - empty cell
  o - omitted
  x - estimated

*/

// After LASSO, we run a flexible parametric model including all these selected variables (by "cross validation")


/*
local predictorvar "age_cat smoking physical_work work_asbestos MetSyn_new ethnic_cat DBP_cat alc_freq_new_2 children_no"
foreach p of local predictorvar {
quietly tab `p', generate(`p') 
}

*/

stset lastfu_prostate, failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    

stpm2 age_cat2 age_cat3 age_cat4 age_cat5 shift_work smoking2 smoking3 physical_work4 physical_work3 physical_work2 ipaq3 ipaq2 ever_PSA crp_cat work_asbestos3 work_asbestos2 father_prostate sibling_prostate wc_cat low_HDL work_per_week_cat marital MetSyn_new5 MetSyn_new4 MetSyn_new3 MetSyn_new2 ethnic_cat3 ethnic_cat2 DBP_cat3 DBP_cat2 fruit_freq firstsex_cat alc_freq_new_25 alc_freq_new_24 alc_freq_new_23 alc_freq_new_22 hba1c_cat work_pesticides children_no3 children_no2 children_no1, scale(h) df(3) eform  


/*

-----------------------------------------------------------------------------------
                  |     exp(b)   Std. Err.      z    P>|z|     [95% Conf. Interval]
------------------+----------------------------------------------------------------
xb                |
         age_cat2 |   3.653531   .7799254     6.07   0.000     2.404392    5.551626
         age_cat3 |   6.379118    1.31236     9.01   0.000     4.262311    9.547203
         age_cat4 |   8.961295   2.075395     9.47   0.000     5.691626    14.10929
         age_cat5 |   12.46375    4.87793     6.45   0.000      5.78781    26.84002
       shift_work |   .6247861   .1207839    -2.43   0.015     .4277369    .9126117
         smoking2 |   .7848517   .1032732    -1.84   0.066     .6064344    1.015761
         smoking3 |   .9486955   .2286238    -0.22   0.827     .5915587    1.521443
   physical_work4 |   1.259457   .4278916     0.68   0.497     .6471336    2.451168
   physical_work3 |   .7423206   .2477216    -0.89   0.372     .3859536    1.427736
   physical_work2 |   1.015613   .1720583     0.09   0.927     .7286574    1.415575
            ipaq3 |   .8248272   .1428011    -1.11   0.266     .5874835    1.158058
            ipaq2 |   1.018684   .1572962     0.12   0.905     .7526695    1.378715
         ever_PSA |   1.277495   .1588441     1.97   0.049     1.001199    1.630038
          crp_cat |   .7370933   .0916291    -2.45   0.014     .5777081    .9404517
   work_asbestos3 |    2.28208    .752448     2.50   0.012     1.195838    4.355011
   work_asbestos2 |   1.242073   .2188731     1.23   0.219     .8793301    1.754456
  father_prostate |   1.875825   .3177022     3.71   0.000     1.345945    2.614313
 sibling_prostate |   2.270626   .7130548     2.61   0.009     1.226984    4.201963
           wc_cat |   .6201315   .1364101    -2.17   0.030     .4029446     .954382
          low_HDL |   .7454551     .17675    -1.24   0.215     .4683803    1.186436
work_per_week_cat |   .9647371   .1631741    -0.21   0.832     .6925292     1.34394
          marital |   .9797368   .3174761    -0.06   0.950     .5191396     1.84899
      MetSyn_new5 |   1.863786   .8519658     1.36   0.173      .760849    4.565556
      MetSyn_new4 |   1.056537   .3322803     0.17   0.861     .5704041    1.956983
      MetSyn_new3 |   1.024216   .2339515     0.10   0.917     .6545762    1.602593
      MetSyn_new2 |   .9039692    .197992    -0.46   0.645     .5884629    1.388635
      ethnic_cat3 |   6.024728   3.565628     3.03   0.002     1.888726     19.2179
      ethnic_cat2 |   1.776414   .6912967     1.48   0.140     .8285077    3.808832
         DBP_cat3 |   1.150756   .1925907     0.84   0.401     .8289455    1.597499
         DBP_cat2 |   1.132722   .1718927     0.82   0.412     .8413015    1.525087
       fruit_freq |   1.178093   .2587718     0.75   0.456      .765968     1.81196
     firstsex_cat |   .9932475   .1404156    -0.05   0.962     .7528751    1.310364
    alc_freq_new5 |    .791242   .2057344    -0.90   0.368     .4753184    1.317146
    alc_freq_new4 |   1.005719   .2449812     0.02   0.981     .6239276    1.621135
    alc_freq_new3 |   .9337404   .1597172    -0.40   0.689     .6677706    1.305645
    alc_freq_new2 |   1.171252   .1733727     1.07   0.286     .8763001    1.565482
        hba1c_cat |   .9653859   .2745332    -0.12   0.901     .5528905    1.685632
  work_pesticides |   1.473399    .409998     1.39   0.164     .8540029    2.542035
     children_no3 |   1.238507    .182394     1.45   0.146     .9279891     1.65293
     children_no2 |   1.444551   .2879302     1.85   0.065     .9773952     2.13499
     children_no1 |   .9565361   .2148813    -0.20   0.843     .6158614    1.485661
            _rcs1 |   1.217334   .0145399    16.47   0.000     1.189167    1.246167
            _rcs2 |   1.040949   .0082864     5.04   0.000     1.024834    1.057317
            _rcs3 |   1.053482   .0062222     8.82   0.000     1.041357    1.065748
            _cons |   .0050399    .002287   -11.66   0.000     .0020709    .0122654
-----------------------------------------------------------------------------------
Note: Estimates are transformed only in the first equation.


*/

// After LASSO, we run a flexible parametric model including all these selected variables (by "adaptive")


stpm2 homosexpartner_no3 homosexpartner_no2 homosex age_cat5 age_cat4 age_cat3 age_cat2 hba1c_cat shift_work DM smoking3 smoking2 ipaq3 ipaq2 MetSyn_new5 MetSyn_new4 MetSyn_new3 MetSyn_new2 SBP_cat4 SBP_cat3 SBP_cat2 wc_cat alc_freq_new_25 alc_freq_new_24 alc_freq_new_23 alc_freq_new_22 physical_work work_asbestos hypertension father_prostate ethnic4 ethnic3 ethnic2 crp_cat DBP_cat3 DBP_cat2 work_per_week_cat work_pesticides hyperlipidemia marital bmi_cat3 bmi_cat2 fruit_freq sibling_prostate, scale(h) df(3) eform  

/*
note: homosex omitted because of collinearity

Iteration 0:   log likelihood =   -1773.83  
Iteration 1:   log likelihood = -1681.5384  
Iteration 2:   log likelihood = -1665.5189  
Iteration 3:   log likelihood =  -1656.551  
Iteration 4:   log likelihood = -1643.9503  
Iteration 5:   log likelihood = -1634.6944  
Iteration 6:   log likelihood = -1634.1343  
Iteration 7:   log likelihood = -1634.1299  
Iteration 8:   log likelihood = -1634.1299  

Log likelihood = -1634.1299                     Number of obs     =     15,836

------------------------------------------------------------------------------------
                   |     exp(b)   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------------+----------------------------------------------------------------
xb                 |
homosexpartner_no3 |   1.397658   .5805672     0.81   0.420     .6191899    3.154844
homosexpartner_no2 |   1.412974   .6414227     0.76   0.446     .5804022    3.439846
           homosex |          1  (omitted)
          age_cat5 |    13.2621   5.081594     6.75   0.000     6.258353    28.10375
          age_cat4 |   10.30012   2.281189    10.53   0.000     6.673058    15.89863
          age_cat3 |   6.993109   1.384349     9.82   0.000     4.744247    10.30797
          age_cat2 |   3.694933   .7782889     6.20   0.000     2.445187    5.583428
         hba1c_cat |   3.912792   4.010261     1.33   0.183     .5249059    29.16702
        shift_work |   .6235094   .1202364    -2.45   0.014     .4272665    .9098863
                DM |    .273056   .2777819    -1.28   0.202     .0371807    2.005328
          smoking3 |   .8053669   .1975703    -0.88   0.378     .4979428    1.302591
          smoking2 |   .7572099   .0984298    -2.14   0.032     .5869059    .9769314
             ipaq3 |   .8007952   .1384257    -1.29   0.199     .5706667    1.123726
             ipaq2 |   1.021301   .1556866     0.14   0.890     .7575247    1.376925
       MetSyn_new5 |   1.143868   .7813096     0.20   0.844     .2998953    4.362969
       MetSyn_new4 |   .7187314   .4104986    -0.58   0.563     .2346435     2.20153
       MetSyn_new3 |   .7850316   .3726613    -0.51   0.610     .3096096    1.990489
       MetSyn_new2 |   .8214808   .2578227    -0.63   0.531     .4440653    1.519665
          SBP_cat4 |   .4912286   .1782304    -1.96   0.050      .241238     1.00028
          SBP_cat3 |   .5635416   .2069807    -1.56   0.118     .2743417    1.157604
          SBP_cat2 |   .6541979   .1568543    -1.77   0.077     .4089028    1.046642
            wc_cat |   .6974806   .1772083    -1.42   0.156     .4239045    1.147615
   alc_freq_new_25 |   1.385108   .3597113     1.25   0.210     .8325823    2.304305
   alc_freq_new_24 |   1.531637   .3925546     1.66   0.096      .926819    2.531142
   alc_freq_new_23 |   1.230515   .3268597     0.78   0.435     .7311112    2.071048
   alc_freq_new_22 |   1.386669   .4323259     1.05   0.294     .7526469    2.554786
     physical_work |    1.01372   .0943975     0.15   0.884     .8446059    1.216695
     work_asbestos |    1.39231    .180949     2.55   0.011     1.079223    1.796225
      hypertension |   1.551194    .621842     1.10   0.273      .707028    3.403262
   father_prostate |   1.929163   .3248804     3.90   0.000     1.386828    2.683586
           ethnic4 |   5.694818   3.368042     2.94   0.003     1.786736    18.15095
           ethnic3 |     1.3076   .7699941     0.46   0.649     .4123213    4.146808
           ethnic2 |   2.736781   1.388884     1.98   0.047     1.012197    7.399713
           crp_cat |   .7271233   .0908775    -2.55   0.011     .5691452    .9289517
          DBP_cat3 |   1.210162   .2379059     0.97   0.332     .8232003    1.779022
          DBP_cat2 |   1.164516   .1923213     0.92   0.356     .8424968    1.609616
 work_per_week_cat |   .9121027   .1553965    -0.54   0.589     .6531652    1.273692
   work_pesticides |   1.435502   .3990574     1.30   0.193     .8324864    2.475316
    hyperlipidemia |   1.165413   .2969768     0.60   0.548     .7072486    1.920382
           marital |   .9376513   .2903481    -0.21   0.835     .5110493    1.720362
          bmi_cat3 |   1.133583   .2807499     0.51   0.613     .6976538    1.841904
          bmi_cat2 |   1.123008   .1612953     0.81   0.419     .8474742    1.488124
        fruit_freq |   1.110931   .2488082     0.47   0.639     .7162243    1.723159
  sibling_prostate |   2.354683   .7379256     2.73   0.006     1.274024    4.351986
             _rcs1 |   1.211582     .01407    16.53   0.000     1.184317    1.239475
             _rcs2 |   1.041817   .0081864     5.21   0.000     1.025895    1.057986
             _rcs3 |   1.053097     .00612     8.90   0.000      1.04117     1.06516
             _cons |   .0054483   .0026082   -10.89   0.000     .0021319    .0139233
------------------------------------------------------------------------------------
Note: Estimates are transformed only in the first equation.
*/



// After LASSO...however we run a flexible parametric model by a traditional methods of adding clinically important variables one by one, including some of the selected variables (by "adaptive") // 
 
 
 
 
  
 
 
****************************// Analysis using classic analysis //************************************

// Codes for flexible parametric survival model (classic analysis)
stset lastfu_prostate, failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    

//univariable in Table 3
local predictorvar "MetSyn_bi DM hypertension hyperlipidemia low_HDL obese "
foreach p of local predictorvar {
stpm2 `p'  if  s_40021_0_0!="SCOT", scale(h) df(4) eform
}
 
// Table 3 model 1 (metabolic syndrome + age + ethnicity + model 1 + deprivation index)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5  if  s_40021_0_0!="SCOT", scale(h) df(4) eform

// Table 3 model 2 (model 1 + smoking + fruit intake + processed meat intake + BMI + exercise level)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 if  s_40021_0_0!="SCOT", scale(h) df(4) eform

// Table 3 model 3 (model 2 + ever had prostate-specific antigen test + father had prostate cancer + sibling had prostate cancer)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate if  s_40021_0_0!="SCOT", scale(h) df(4) eform

// Table 3 model 4 (model 3 + C-reactive protein level + Insulin-like growth factor and testosterone)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone  if  s_40021_0_0!="SCOT", scale(h) df(4) eform
 
// Individual metabolic syndrome conditions - can put in Supplementary Figure
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone  if  s_40021_0_0!="SCOT", scale(h) df(4) eform


 
 
****************************// Analysis using IPTW //************************************

// Codes for Estimating the Inverse Probability of Treatment Weights (IPTW)   
logistic MetSyn_bi ipaq age i.ethnic bmi medications fruit_freq alc_freq_new_2 processedmeat_freq i.smoking, nolog       
predict double ps             //ps = propensity score
gen double HAW = ((MetSyn_bi == 1) / ps) + ((MetSyn_bi == 0) / (1 - ps))     // Compute the inverse probability Treatment weights (IPTW)
summarize HAW, detail

stset lastfu_prostate [pw=HAW], failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    

//univariable in Table 3
local predictorvar "MetSyn_bi DM hypertension hyperlipidemia low_HDL obese i.age_cat marital i.ethnic i.qdi i.smoking2 i.processedmeat_freq i.fruit_freq i.ipaq bmi i.ever_PSA i.father_prostate i.sibling_prostate i.crp_cat testosterone IGF"
foreach p of local predictorvar {
stpm2 `p'  if  s_40021_0_0!="SCOT", scale(h) df(4) eform
}



// Table 3 model 1 (metabolic syndrome + age + ethnicity + model 1 + deprivation index)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5  if  s_40021_0_0!="SCOT", scale(h) df(4) eform

// Table 3 model 2 (model 1 + smoking + fruit intake + processed meat intake + BMI + exercise level)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 if  s_40021_0_0!="SCOT", scale(h) df(4) eform

// Table 3 model 3 (model 2 + ever had prostate-specific antigen test + father had prostate cancer + sibling had prostate cancer)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate if  s_40021_0_0!="SCOT", scale(h) df(4) eform

// Table 3 model 4 (model 3 + C-reactive protein level + Insulin-like growth factor and testosterone)
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone  if  s_40021_0_0!="SCOT", scale(h) df(4) eform
 
// Individual metabolic syndrome conditions - can put in Supplementary Figure
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone  if  s_40021_0_0!="SCOT", scale(h) df(4) eform

 
// Below are Standsurv codes
stset competing_risk_date_new [pw=HAW], failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store prostate_incidence

stset competing_risk_date_new [pw=HAW], failure(prostate==0) origin(time assessment_date) id(id) scale(365.25)    
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5  ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store comp_incidence

range tt 0 15 30

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence comp_incidence) cif ci timevar(tt) atvar(F) verbose 
// Plot prostate_incidence 
twoway  (rarea F_prostate_incidence_lci F_prostate_incidence_uci tt, color(red%30)) ///
        (line F_prostate_incidence tt, color(red)) ///
		(rarea F_comp_incidence_lci F_comp_incidence_uci tt, color(blue%30)) ///
        (line F_comp_incidence tt, color(blue)) ///
        , legend(order(2 "Prostate Incidence"  4 "No prostate cancer") cols(1) ring(0) pos(11)) ///
		ylabel(,angle(h) format(%3.2f)) ///
		xtitle("Follow-up (years)") ytitle("Cumulative incidence") ///
		name(cifs, replace) ///
		saving(prostate_incidence, replace) nodraw
			
	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence comp_incidence) cif ci timevar(tt) verbose ///
    at1(MetSyn_bi 1) at2(MetSyn_bi 0) atvars(F_Met1 F_Met2)  atref(2) ///
    contrast(ratio) contrastvars(r1)   //ratio
	

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence comp_incidence) cif ci timevar(tt) verbose ///
    at1(MetSyn_bi 1) at2(MetSyn_bi 0) atvars(F_Met1_d F_Met2_d)  atref(2) ///
    contrast(difference) contrastvars(d1)  //difference
	
// Plot prostate cancer incidence based on metabolic syndrome
twoway  (rarea F_Met1_d_prostate_incidence_lci F_Met1_d_prostate_incidence_uci tt, color(red%30)) ///
        (line F_Met1_d_prostate_incidence tt, color(red)) ///
        (rarea F_Met2_d_prostate_incidence_lci F_Met2_d_prostate_incidence_uci tt, color(blue%30)) ///
		(line F_Met2_d_prostate_incidence tt, color(blue)) ///
         , legend(order(2 "Metabolic Syndrome" 4 "No Metabolic Syndrome") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-up (years)") ytitle("Cumulative incidence") ///
		name(metabolic, replace) ///
		saving(metabolic, replace) nodraw

// Plot on competing risks based on metabolic syndrome
twoway  (rarea F_Met1_comp_incidence_lci F_Met1_comp_incidence_uci tt, color(red%30)) ///
        (line F_Met1_comp_incidence tt, color(red)) ///
		(rarea F_Met2_comp_incidence_lci F_Met2_comp_incidence_uci tt, color(blue%30)) ///
        (line F_Met2_comp_incidence tt, color(blue)) ///
        , legend(order(1 "Metabolic Syndrome" 3 "No Metabolic Syndrome") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic, replace) ///
		saving(metabolic_othercancer_death, replace) nodraw

stset competing_risk_date_new [pw=HAW], failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store prostate_incidence2

stset competing_risk_date_new [pw=HAW], failure(prostate==0) origin(time assessment_date) id(id) scale(365.25)    
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5  ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store comp_incidence2


// Plot prostate cancer incidence based on individual metabolic syndrome component - DM
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose at1(DM 1) at2(DM 0) atvars(F_DM1 F_DM2) atref(2) contrast(ratio) contrastvars(R_DM)

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose at1(DM 1) at2(DM 0) atvars(F_DM1 F_DM2) atref(2) contrast(difference) contrastvars(D_DM)
	
twoway  (rarea F_DM1_prostate_incidence2_lci F_DM1_prostate_incidence2_uci tt, color(red%30)) ///
        (line F_DM1_prostate_incidence2 tt, color(red)) ///
        (rarea F_DM2_prostate_incidence2_lci F_DM2_prostate_incidence2_uci tt, color(blue%30)) ///
		(line F_DM2_prostate_incidence2 tt, color(blue)) ///
         , legend(order(2 "Diabetes" 4 "No Diabetes") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic_DM, replace) ///
		saving(metabolic_DM, replace) nodraw
	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(obese 1) at2(obese 0) atvars(F_ob1 F_ob2)  atref(2) ///
    contrast(ratio) contrastvars(R_ob)

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(obese 1) at2(obese 0) atvars(F_ob1 F_ob2)  atref(2) ///
    contrast(difference) contrastvars(D_ob)
	
// Plot prostate cancer incidence based on individual metabolic syndrome component - obesity
twoway  (rarea F_ob1_prostate_incidence2_lci F_ob1_prostate_incidence2_uci tt, color(red%30)) ///
        (line F_ob1_prostate_incidence2 tt, color(red)) ///
        (rarea F_ob2_prostate_incidence2_lci F_ob2_prostate_incidence2_uci tt, color(blue%30)) ///
		(line F_ob2_prostate_incidence2 tt, color(blue)) ///
         , legend(order(2 "Obesity" 4 "No Obesity") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic_obesity, replace) ///
		saving(metabolic_obesity, replace) nodraw
						
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(hypertension 1) at2(hypertension 0) atvars(F_HT1 F_HT2)  atref(2) ///
    contrast(ratio) contrastvars(r_HT)

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(hypertension 1) at2(hypertension 0) atvars(F_HT1 F_HT2)  atref(2) ///
    contrast(difference) contrastvars(D_HT)
	
// Plot prostate cancer incidence based on individual metabolic syndrome component - hypertension
twoway  (rarea F_HT1_prostate_incidence2_lci F_HT1_prostate_incidence2_uci tt, color(red%30)) ///
        (line F_HT1_prostate_incidence2 tt, color(red)) ///
        (rarea F_HT2_prostate_incidence2_lci F_HT2_prostate_incidence2_uci tt, color(blue%30)) ///
		(line F_HT2_prostate_incidence2 tt, color(blue)) ///
         , legend(order(2 "Hypertension" 4 "No Hypertension") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic_hypertension, replace) ///
		saving(metabolic_hypertension, replace) nodraw

				
	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(hyperlipidemia 1) at2(hyperlipidemia 0) atvars(F_HL1 F_HL2)  atref(2) ///
    contrast(ratio) contrastvars(r_HL)
	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(hyperlipidemia 1) at2(hyperlipidemia 0) atvars(F_HL1 F_HL2)  atref(2) ///
    contrast(difference) contrastvars(D_HL)
	
// Plot prostate cancer incidence based on individual metabolic syndrome component - hyperlipidemia
twoway  (rarea F_HL1_prostate_incidence2_lci F_HL1_prostate_incidence2_uci tt, color(red%30)) ///
        (line F_HL1_prostate_incidence2 tt, color(red)) ///
        (rarea F_HL2_prostate_incidence2_lci F_HL2_prostate_incidence2_uci tt, color(blue%30)) ///
		(line F_HL2_prostate_incidence2 tt, color(blue)) ///
         , legend(order(2 "Hyperlipidemia" 4 "No Hyperlipidemia") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic_hyperlipidemia, replace) ///
		saving(metabolic_hyperlipidemia, replace) nodraw

			
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(low_HDL 1) at2(low_HDL 0) atvars(F_low1 F_low2)  atref(2) ///
    contrast(ratio) contrastvars(r_low)
	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(low_HDL 1) at2(low_HDL 0) atvars(F_low1 F_low2)  atref(2) ///
    contrast(difference) contrastvars(D_low)
	
	
// Plot prostate cancer incidence based on individual metabolic syndrome component - low HDL level
twoway  (rarea F_low1_prostate_incidence2_lci F_low1_prostate_incidence2_uci tt, color(red%30)) ///
        (line F_low1_prostate_incidence2 tt, color(red)) ///
        (rarea F_low2_prostate_incidence2_lci F_low2_prostate_incidence2_uci tt, color(blue%30)) ///
		(line F_low2_prostate_incidence2 tt, color(blue)) ///
         , legend(order(2 "Low HDL" 4 "HDL Not Low") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic_lowHDL, replace) ///
		saving(metabolic_lowHDL, replace) nodraw
restore






// Below are Standsurv codes
stset competing_risk_otherdeath_date [pw=HAW], failure(competing_risk_otherdeath==1) origin(time assessment_date) id(id) scale(365.25)    
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store prostate_incidence

stset competing_risk_otherdeath_date [pw=HAW], failure(competing_risk_otherdeath==2) origin(time assessment_date) id(id) scale(365.25)    
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5  ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store comp_incidence

stset competing_risk_otherdeath_date [pw=HAW], failure(competing_risk_otherdeath==1) origin(time assessment_date) id(id) scale(365.25)    
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5 marital ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store prostate_incidence2

stset competing_risk_otherdeath_date [pw=HAW], failure(competing_risk_otherdeath==2) origin(time assessment_date) id(id) scale(365.25)    
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5  ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store comp_incidence2


standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & marital!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence comp_incidence) cif ci timevar(tt) verbose ///
    at1(MetSyn_bi 1) at2(MetSyn_bi 0) atvars(F_Met1_d2 F_Met2_d2)  atref(2) ///
    contrast(difference) contrastvars(d12)  //difference


standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & marital!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose at1(DM 1) at2(DM 0) atvars(F_DM12 F_DM22) atref(2) contrast(difference) contrastvars(D_DM2)


standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & marital!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(obese 1) at2(obese 0) atvars(F_ob12 F_ob22)  atref(2) ///
    contrast(difference) contrastvars(D_ob2)

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & marital!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(hypertension 1) at2(hypertension 0) atvars(F_HT12 F_HT22)  atref(2) ///
    contrast(difference) contrastvars(D_HT2)	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & marital!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(hyperlipidemia 1) at2(hyperlipidemia 0) atvars(F_HL12 F_HL22)  atref(2) ///
    contrast(difference) contrastvars(D_HL2)
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & marital!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence2 comp_incidence2) cif ci timevar(tt) verbose ///
    at1(low_HDL 1) at2(low_HDL 0) atvars(F_low12 F_low22)  atref(2) ///
    contrast(difference) contrastvars(D_low2)

list F_Met1_prostate_incidence F_Met1_prostate_incidence_lci F_Met1_prostate_incidence_uci if tt==5, noobs ab(30)  // metabolic syndrome
list F_Met2_prostate_incidence F_Met2_prostate_incidence_lci F_Met2_prostate_incidence_uci if tt==1, noobs ab(30)  // no metabolic syndrome
list r1_prostate_incidence r1_prostate_incidence_lci r1_prostate_incidence_uci tt if tt<=5 , noobs ab(30)   //ratio of metabolic syndrome vs no metabolic syndrome 5 years
list r1_prostate_incidence r1_prostate_incidence_lci r1_prostate_incidence_uci tt if tt<=10 , noobs ab(30)   //ratio of metabolic syndrome vs no metabolic syndrome 10 years
list r1_prostate_incidence r1_prostate_incidence_lci r1_prostate_incidence_uci tt if tt<=15 , noobs ab(30)   //ratio of metabolic syndrome vs no metabolic syndrome 15 years



//Previous analysis//
list d1_prostate_incidence d1_prostate_incidence_lci d1_prostate_incidence_uci tt if tt<=1 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 1 years
list d1_prostate_incidence d1_prostate_incidence_lci d1_prostate_incidence_uci tt if tt<=3 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 3 years
list d1_prostate_incidence d1_prostate_incidence_lci d1_prostate_incidence_uci tt if tt<=5 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 5 years
list d1_prostate_incidence d1_prostate_incidence_lci d1_prostate_incidence_uci tt if tt<=10 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 10 years

list D_DM_prostate_incidence2 D_DM_prostate_incidence2_lci D_DM_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of DM vs no DM 1 years
list D_DM_prostate_incidence2 D_DM_prostate_incidence2_lci D_DM_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of DM vs no DM 3 years
list D_DM_prostate_incidence2 D_DM_prostate_incidence2_lci D_DM_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of DM vs no DM 5 years
list D_DM_prostate_incidence2 D_DM_prostate_incidence2_lci D_DM_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of DM vs no DM 10 years

list D_ob_prostate_incidence2 D_ob_prostate_incidence2_lci D_ob_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of obesity vs no obesity 1 years
list D_ob_prostate_incidence2 D_ob_prostate_incidence2_lci D_ob_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of obesity vs no obesity 3 years
list D_ob_prostate_incidence2 D_ob_prostate_incidence2_lci D_ob_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of obesity vs no obesity 5 years
list D_ob_prostate_incidence2 D_ob_prostate_incidence2_lci D_ob_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of obesity vs no obesity 10 years

list D_HT_prostate_incidence2 D_HT_prostate_incidence2_lci D_HT_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of hypertension vs no hypertension 1 years
list D_HT_prostate_incidence2 D_HT_prostate_incidence2_lci D_HT_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of hypertension vs no hypertension 3 years
list D_HT_prostate_incidence2 D_HT_prostate_incidence2_lci D_HT_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of hypertension vs no hypertension 5 years
list D_HT_prostate_incidence2 D_HT_prostate_incidence2_lci D_HT_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of hypertension vs no hypertension 10 years

list D_HL_prostate_incidence2 D_HL_prostate_incidence2_lci D_HL_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of hyperlipidemia vs no hyperlipidemia 1 years
list D_HL_prostate_incidence2 D_HL_prostate_incidence2_lci D_HL_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of hyperlipidemia vs no hyperlipidemia 3 years
list D_HL_prostate_incidence2 D_HL_prostate_incidence2_lci D_HL_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of hyperlipidemia vs no hyperlipidemia 5 years
list D_HL_prostate_incidence2 D_HL_prostate_incidence2_lci D_HL_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of hyperlipidemia vs no hyperlipidemia 10 years

list D_low_prostate_incidence2 D_low_prostate_incidence2_lci D_low_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of low HDL-C vs no low HDL-C 1 years
list D_low_prostate_incidence2 D_low_prostate_incidence2_lci D_low_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of low HDL-C vs no low HDL-C 3 years
list D_low_prostate_incidence2 D_low_prostate_incidence2_lci D_low_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of low HDL-C vs no low HDL-C 5 years
list D_low_prostate_incidence2 D_low_prostate_incidence2_lci D_low_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of low HDL-C vs no low HDL-C 10 years


// Analysis on 16/8/2023//
list d12_prostate_incidence d12_prostate_incidence_lci d12_prostate_incidence_uci tt if tt<=1 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 1 years
list d12_prostate_incidence d12_prostate_incidence_lci d12_prostate_incidence_uci tt if tt<=3 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 3 years
list d12_prostate_incidence d12_prostate_incidence_lci d12_prostate_incidence_uci tt if tt<=5 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 5 years
list d12_prostate_incidence d12_prostate_incidence_lci d12_prostate_incidence_uci tt if tt<=10 , noobs ab(30)   //Risk difference of metabolic syndrome vs no metabolic syndrome 10 years

list D_DM2_prostate_incidence2 D_DM2_prostate_incidence2_lci D_DM2_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of DM vs no DM 1 years
list D_DM2_prostate_incidence2 D_DM2_prostate_incidence2_lci D_DM2_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of DM vs no DM 3 years
list D_DM2_prostate_incidence2 D_DM2_prostate_incidence2_lci D_DM2_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of DM vs no DM 5 years
list D_DM2_prostate_incidence2 D_DM2_prostate_incidence2_lci D_DM2_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of DM vs no DM 10 years

list D_ob2_prostate_incidence2 D_ob2_prostate_incidence2_lci D_ob2_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of obesity vs no obesity 1 years
list D_ob2_prostate_incidence2 D_ob2_prostate_incidence2_lci D_ob2_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of obesity vs no obesity 3 years
list D_ob2_prostate_incidence2 D_ob2_prostate_incidence2_lci D_ob2_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of obesity vs no obesity 5 years
list D_ob2_prostate_incidence2 D_ob2_prostate_incidence2_lci D_ob2_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of obesity vs no obesity 10 years

list D_HT2_prostate_incidence2 D_HT2_prostate_incidence2_lci D_HT2_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of hypertension vs no hypertension 1 years
list D_HT2_prostate_incidence2 D_HT2_prostate_incidence2_lci D_HT2_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of hypertension vs no hypertension 3 years
list D_HT2_prostate_incidence2 D_HT2_prostate_incidence2_lci D_HT2_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of hypertension vs no hypertension 5 years
list D_HT2_prostate_incidence2 D_HT2_prostate_incidence2_lci D_HT2_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of hypertension vs no hypertension 10 years

list D_HL2_prostate_incidence2 D_HL2_prostate_incidence2_lci D_HL2_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of hyperlipidemia vs no hyperlipidemia 1 years
list D_HL2_prostate_incidence2 D_HL2_prostate_incidence2_lci D_HL2_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of hyperlipidemia vs no hyperlipidemia 3 years
list D_HL2_prostate_incidence2 D_HL2_prostate_incidence2_lci D_HL2_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of hyperlipidemia vs no hyperlipidemia 5 years
list D_HL2_prostate_incidence2 D_HL2_prostate_incidence2_lci D_HL2_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of hyperlipidemia vs no hyperlipidemia 10 years

list D_low2_prostate_incidence2 D_low2_prostate_incidence2_lci D_low2_prostate_incidence2_uci tt if tt<=1 , noobs ab(30)   //Risk difference of low HDL-C vs no low HDL-C 1 years
list D_low2_prostate_incidence2 D_low2_prostate_incidence2_lci D_low2_prostate_incidence2_uci tt if tt<=3 , noobs ab(30)   //Risk difference of low HDL-C vs no low HDL-C 3 years
list D_low2_prostate_incidence2 D_low2_prostate_incidence2_lci D_low2_prostate_incidence2_uci tt if tt<=5 , noobs ab(30)  //Risk difference of low HDL-C vs no low HDL-C 5 years
list D_low2_prostate_incidence2 D_low2_prostate_incidence2_lci D_low2_prostate_incidence2_uci tt if tt<=10 , noobs ab(30)   //Risk difference of low HDL-C vs no low HDL-C 10 years








//**** Death due to causes other than prostate cancer*******

gen competing_risk_otherdeath = 0    // 0 = no prostate cancer or other cancer, no death; 1 = prostate cancer; 2 = other cancers or any death
replace competing_risk_otherdeath = 1 if prostate ==1 
replace competing_risk_otherdeath = 2 if death==1 & prostate ==0

gen competing_risk_otherdeath_date = date("2019, 7, 31", "YMD")
replace competing_risk_otherdeath_date = DOD if competing_risk_otherdeath==2
replace competing_risk_otherdeath_date = prostate_date if competing_risk_otherdeath==1
format %tdDD/NN/CCYY competing_risk_otherdeath_date

tab competing_risk_otherdeath if s_40021_0_0!="SCOT" & sex==1



// Below are Standsurv codes to plot the "causes other than prostate cancer"
stset competing_risk_otherdeath_date [pw=HAW], failure(competing_risk_otherdeath==1) origin(time assessment_date) id(id) scale(365.25)    
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store prostate_incidence

stset competing_risk_otherdeath_date [pw=HAW], failure(competing_risk_otherdeath==2) origin(time assessment_date) id(id) scale(365.25)    
stpm2 MetSyn_bi age_cat2 age_cat3 age_cat4 age_cat5 ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store comp_incidence

range tt 0 15 30

standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence comp_incidence) cif ci timevar(tt) atvar(F) verbose 
// Plot prostate_incidence 
twoway  (rarea F_prostate_incidence_lci F_prostate_incidence_uci tt, color(red%30)) ///
        (line F_prostate_incidence tt, color(red)) ///
		(rarea F_comp_incidence_lci F_comp_incidence_uci tt, color(blue%30)) ///
        (line F_comp_incidence tt, color(blue)) ///
        , legend(order(2 "Prostate Incidence"  4 "Other cancers or death") cols(1) ring(0) pos(11)) ///
		ylabel(,angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(cifs, replace) ///
		saving(prostate_incidence, replace) nodraw
			
	
standsurv if ethnic4!=. & smoking3!=. & qdi2!=. & fruit_freq!=. & ever_PSA!=. & father_prostate!=. & processedmeat_freq3!=. & sibling_prostate!=. & ipaq3!=. & crp_cat!=. & bmi!=. & IGF!=. & testosterone!=., crmodels(prostate_incidence comp_incidence) cif ci timevar(tt) verbose ///
    at1(MetSyn_bi 1) at2(MetSyn_bi 0) atvars(F_Met1 F_Met2)  atref(2) ///
    contrast(ratio) contrastvars(r1)
	
// Plot prostate cancer incidence based on metabolic syndrome
twoway  (rarea F_Met1_prostate_incidence_lci F_Met1_prostate_incidence_uci tt, color(red%30)) ///
        (line F_Met1_prostate_incidence tt, color(red)) ///
        (rarea F_Met2_prostate_incidence_lci F_Met2_prostate_incidence_uci tt, color(blue%30)) ///
		(line F_Met2_prostate_incidence tt, color(blue)) ///
         , legend(order(2 "Metabolic Syndrome" 4 "No Metabolic Syndrome") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic, replace) ///
		saving(metabolic, replace) nodraw

// Plot on competing risks based on metabolic syndrome
twoway  (rarea F_Met1_comp_incidence_lci F_Met1_comp_incidence_uci tt, color(red%30)) ///
        (line F_Met1_comp_incidence tt, color(red)) ///
		(rarea F_Met2_comp_incidence_lci F_Met2_comp_incidence_uci tt, color(blue%30)) ///
        (line F_Met2_comp_incidence tt, color(blue)) ///
        , legend(order(1 "Metabolic Syndrome" 3 "No Metabolic Syndrome") cols(1) ring(0) pos(11)) ///
		ylabel(, angle(h) format(%3.2f)) ///
		xtitle("Follow-Up (Years)") ytitle("Adjusted cause-specific cumulative incidence") ///
		name(metabolic, replace) ///
		saving(metabolic_othercancer_death, replace) nodraw


stset competing_risk_date_new [pw=HAW], failure(prostate==1) origin(time assessment_date) id(id) scale(365.25)    
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5 ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store prostate_incidence2

stset competing_risk_date_new [pw=HAW], failure(prostate==0) origin(time assessment_date) id(id) scale(365.25)    
stpm2 obese DM hypertension hyperlipidemia low_HDL age_cat2 age_cat3 age_cat4 age_cat5 ethnic4 ethnic3 ethnic2 qdi2 qdi3 qdi4 qdi5 smoking3 smoking2 fruit_freq processedmeat_freq3 processedmeat_freq2 bmi ipaq3 ipaq2 ever_PSA father_prostate sibling_prostate crp_cat IGF testosterone if  s_40021_0_0!="SCOT", scale(h) df(4) eform
estimates store comp_incidence2
























/* Assessing IPTW overlap by hand: metabolic syndrome*/
sort MetSyn_bi
by MetSyn_bi: summarize ps
// propsensity score
// Same for the ipcw weights
kdensity ps if MetSyn_bi==1, generate(x1pointsa d1A) nograph n(10000)
kdensity ps if MetSyn_bi==0, generate(x0pointsa d0A) nograph n(10000)
label variable d1A "density for metabolic syndrome=1"
label variable d0A "density for metabolic syndrome=0"
twoway (line d0A x0pointsa , yaxis(1))(line d1A x1pointsa, yaxis(2)), name(iptw_metabolic, replace) saving(iptw_metabolic, replace) nodraw
    
/* Assessing IPTW overlap by hand: Obesity*/
sort obese
by obese: summarize ps
// propsensity score
// Same for the ipcw weights
kdensity ps if obese==1, generate(x1pointsa_obese d1A_obese) nograph n(10000)
kdensity ps if obese==0, generate(x0pointsa_obese d0A_obese) nograph n(10000)
label variable d1A_obese "density for obesity=1"
label variable d0A_obese "density for obesity=0"
twoway (line d0A_obese x0pointsa_obese, yaxis(1))(line d1A_obese x1pointsa_obese, yaxis(2)), name(iptw_obese, replace) saving(iptw_obese, replace) nodraw
    
/* Assessing IPTW overlap by hand: Diabetes*/
sort DM
by DM: summarize ps
// propsensity score
// Same for the ipcw weights
kdensity ps if DM==1, generate(x1pointsa_DM d1A_DM) nograph n(10000)
kdensity ps if DM==0, generate(x0pointsa_DM d0A_DM) nograph n(10000)
label variable d1A_DM "density for diabetes=1"
label variable d0A_DM "density for diabetes=0"
twoway (line d0A_DM x0pointsa_DM, yaxis(1))(line d1A_DM x1pointsa_DM, yaxis(2)), name(iptw_DM, replace) saving(iptw_DM, replace) nodraw
   
/* Assessing IPTW overlap by hand: Hypertension*/
sort hypertension
by hypertension: summarize ps
// propsensity score
// Same for the ipcw weights
kdensity ps if hypertension==1, generate(x1pointsa_HT d1A_HT) nograph n(10000)
kdensity ps if hypertension==0, generate(x0pointsa_HT d0A_HT) nograph n(10000)
label variable d1A_HT "density for hypertension=1"
label variable d0A_HT "density for hypertension=0"
twoway (line d0A_HT x0pointsa_HT, yaxis(1))(line d1A_HT x1pointsa_HT, yaxis(2)), name(iptw_HT, replace) saving(iptw_HT, replace) nodraw
    
/* Assessing IPTW overlap by hand: Hyperlipidemia*/
sort hyperlipidemia
by hyperlipidemia: summarize ps
// propsensity score
// Same for the ipcw weights
kdensity ps if hyperlipidemia==1, generate(x1pointsa_HL d1A_HL) nograph n(10000)
kdensity ps if hyperlipidemia==0, generate(x0pointsa_HL d0A_HL) nograph n(10000)
label variable d1A_HL "density for Hyperlipidemia=1"
label variable d0A_HL "density for Hyperlipidemia=0"
twoway (line d0A_HL x0pointsa_HL , yaxis(1))(line d1A_HL x1pointsa_HL, yaxis(2)), name(iptw_HL, replace) saving(iptw_HL, replace) nodraw
    
/* Assessing IPTW overlap by hand: low_HDL*/
sort low_HDL
by low_HDL: summarize ps
// propsensity score
// Same for the ipcw weights
kdensity ps if low_HDL==1, generate(x1pointsa_lowHDL d1A_lowHDL) nograph n(10000)
kdensity ps if low_HDL==0, generate(x0pointsa_lowHDL d0A_lowHDL) nograph n(10000)
label variable d1A_lowHDL "density for low HDL level=1"
label variable d0A_lowHDL "density for low HDL level=0"
twoway (line d0A_lowHDL x0pointsa_lowHDL, yaxis(1))(line d1A_lowHDL x1pointsa_lowHDL, yaxis(2)), name(iptw_lowHDL, replace) saving(iptw_lowHDL, replace) nodraw
    
