# Generalities/multiple cancer types
1. _Tumour cell proliferation kinetics and tumour growth rate_ (https://www.tandfonline.com/doi/abs/10.3109/02841868909111193)
    - Doubling times (DT) in human tumours seemed constant in the clinic, except for when longer-term follow-up was possible => there *is* an extension of DT (logistic or Gompertzian growth, not exponential) but this is not picked up when monitoring times are typically around a few doublings.
    - "In summary despite the occurrence of irregular growth patterns, the growth rate is in most patients either constant (exponential growth) or progressively decreasing."
    - **"Frequency distribution of individual values for the DT is log-normal."**
    - Table 1 shows doubling times by histological type-embryonal tumour (27 days), malignant lymphoma (29 days), mesenchymal sarcoma (41 days), squamous cell carcinoma (58 days) and adenocarcinoma (83 days)-in patients with pulmonary metastases; doubling times are log-normal within a histological type and less differentiated types have faster growth.
2. _The temporal order of genetic pathway alterations in tumourigenesis_ (http://cancerpreventionresearch.aacrjournals.org/lookup/doi/10.1158/1940-6207.CAPR-10-0374)
    - Like the Sprouffske et al. paper under _Colon_ but using actual mutation data
    - Demonstrates that order constraints are stronger at the pathway level than the genetic level; potential orders for pathway changes predicted for colorectal, pancreatic and primary glioblastoma


# Skin
1. _Stem cell competition orchestrates skin homeostasis and ageing_ (https://www.nature.com/articles/s41586-019-1085-7)
    - Experimental demonstration of competition within the stem cell pool between high and low expression of the desmosome component gene, _Col17a1_. _In vivo_ clonal analysis in mice and _in vitro_ 3D cultures of human keratinocytes
    - High-expressing stem cells divide symmetrically and spread over the basement membrane while low-expressing stem cells divide asymmetrically. The low expression also impairs their adhesion to the basement membrane and surrounding cells. This leads to delamination of the low-expressing cells and their replacement by high-expressing ones.
2. _Distinct modes of cell competition shape mammalian tissue morphogenesis_ (https://www.nature.com/articles/s41586-019-1199-y)
    - Experimental study of cell competition in mouse embryos
    - _Myc_ is known to play a role in competitive elimination of "loser" cells in Dros imaginal discs-attempt to see how this pans out in mammalian systems using various genetic tools and live imaging
    - _Myc_ does seem to play early in mouse embryo development-_Mycn_-floxed and knockout cells did worse than wildtype cells, grew slower and dead more but only in the neighbourhood of the wildtype in early embryos when the basal layer is being established => some sort of cell-cell contact is necesary. _Mycn_ loss is also not cell autonomous, and win/lose fate depends on the genotypes of the neighbours-when _Mycn_ is knocked out completely, the epidermis is intact, while mosaic knockouts which leave some wildtypes lead to death of the _Mycn_ knockouts.
    - If competition-induced apoptosis is blocked,"loser" cells may be cleared by differentiation-in late-stage skin when strata have been established, more of the loser cells are found in the upper layers than around the basal layer.
    - Physiological consequences for barrier function

# Breast
1. _A Gompertzian model of human breast cancer growth_ (https://cancerres.aacrjournals.org/sites/all/libraries/pdfjs/web/viewer.html?file=/content/canres/48/24_Part_1/7067.full.pdf)
    - Gompertzian model fits to three different breast cancer datasets, showing that with some adjustments and assumptions, it is sufficiently most of the patterns including mortality, time from initiation to clinical detection and time to relapse following surgical intervention
    - Fits this equation: N(t) = N(0).exp[k(1-exp(-bt))]; this separates the growth rate from the asymptote, but probably not the initial density and the asymptote or the initial density from the initial growth rate. This equation is the exact same as the Zweifel and Lasker re-parameterisation of a W_0 form Gompertzian.
    - Growth curves are cell number vs time in years.
    - *b* could be taken as a measure of growth rate, and from the data, they find it to be lognormally-distributed with mean = -2.9 and SD = 0.71.
2. _Quantitation and Gompertzian analysis of tumour growth_ (https://link.springer.com/article/10.1023/A:1005906900231)
    - Tracking of human breast cancer xenografts in nude mice
    - Data of tumour volume against days since implantation are fitted with a transformed Gompertz equation-this thing is a bit weird, seems to combine an asymptote parameter with a growth rate and a separate shape parameter.
    - No values for growth rate given, could be inferred with some effort-reliability uncertain.
3. _Differences in predictions of ODE models of tumour growth: a cautionary example_ (https://bmccancer.biomedcentral.com/articles/10.1186/s12885-016-2164-x) **and** _Best fitting tumour growth models of the von Bertalanffy-Putter type_ (https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5911-y)
    - Two related papers using the same data of tumour volume vs time in days from a metastatic human breast cancer xenograft in nude mice
    - The first paper would wide differences in the predictions of the seven different ODE models from exponential to von Bertalanffy in terms of timelines of growth and regression under therapy.
    - The second paper found the same differences, but also shows improvement in model predictions as more of the longitudinal data is added to the model; this is obvious though-more data, better predictions.
    - The general equation used in the second paper is of the form: dv(t)/dt = p.[v(t)^a] - q.[v(t)^b]
    - The claim is that for different values of a and b, various standard models can be recovered like the logistic (a=1, b=2) and exponential (a=1, b=0); Gompertz requires a different functional form and therefore not included.
    - Table 2 in the second paper gives optimal value of the exponents and the free parameters p, q, and v found by their custom-made "simulated annealing" optimisation process. Utility uncertain.
4. _An experimental-mathematical approach to predict tumor cell growth as a function of glucose availability in breast cancer cell lines_ (https://www.biorxiv.org/content/biorxiv/early/2020/10/05/2020.10.05.326041.full.pdf?%3Fcollection=)
    - *In vitro* measurement of proliferation, death rate and glucose uptake rates for two different breast cancer cell lines-one of HER2+ and the other, triple-negative
    - Proliferation and death rates both in per day, while glucose uptake rates are in mM per cell per day.
5. _GROWTH RATE, KINETICS OF TUMOR CELL PROLIFERATION AND LONG-TERM OUTCOME IN HUMAN BREAST CANCER_ (https://cyber.sci-hub.st/MTAuMTAwMi9pamMuMjkxMDQ0MDEwNA==/tubiana1989.pdf#view=FitH)
    - Labelling Index from tritiated thymidine assays used as a proxy for proliferation rate in primary clinical samples from mammograms
    - LI has a range from < 0.25 to > 3.84, while media tumour volume doubling time was found to be 7 months.
    - Patients with higher LI or lower doubling time had a higher risk of metastasis.

# Prostate
1. _Mathematical modeling of prostate cancer progression in response to androgen ablation therapy_ (https://www.pnas.org/content/108/49/19701)
    - Three cell types, as in Cunningham-normal epithelium (E), androgen-dependent cancer cells (N) and castration-resistant cancer cells (M).
    - N are assumed to have higher net doubling rate, but the value is estimated from patient data, and M have different turnover rates from N. Doubling time and mutation frequency in M is used to explore levels of aggressiveness.
    - Fits of cycles of PSA from patient data; first cycle leads to doubling time of 80 days for N, and second cycle leads to a range of 75-135 days for M cells.
    - Competition between N and M is modelled with a single parameter, theta; for effect of N on M, the term becomes theta*N and for the other, M(1\theta). Theta always > 1.
    - SI Text-
        - Proliferation and apoptosis rates of normal prostate epithelium is from experimental data in SD rats; homeostatic proliferation rate due to androgen binding is calculated as ~ 0.029 per hour or a doubling time of about 35 days. It's curious that this is much higher than the estimates for N and M. Death rate is taken to be constant at 0.0083, for which the source is a different paper listed next. This entire exercise seems to be trying to decompose the measured net proliferation rates into androgen-drive terms that are multiplied by some intrinsic rate of growth. It's this latter rate that is calculated to be 0.029.
        - Testosterone-driven vs DHT-driven effects are different on proliferation and death rates; the former affects only proliferation rates while the latter increases proliferation rate **and** decreases death rate due to some pro-survival signaling implications of DHT.
        - To go from rat to human prostate, the size and time to maturation differences lead to a much higher number of epithelial cells in the mature human prostate. Most kinetic parameters are kept the same. Given that they have an estimate of the "cell-intrinsic" growth rate from the rat data, they can probably calculate the actual proliferation rate from the same equation, but this is speculation-the final value for human prostate epithelial cell proliferation rate is not stated explicitly.
2. _Implications of cell kinetic changes during the progression of human prostatic cancer_ (https://clincancerres.aacrjournals.org/content/clincanres/1/5/473.full.pdf)
    - Various clinical samples from patients undergoing treatment for prostate cancer, with normal tissue, precancer neoplastic cells and cancer cells, as well as bone and soft tissue metastases; direct staining of tissue sections and primary cell cultures _in vitro_
    - Time between successive mitoses found by observation of primary cultures is not very different from that of established prostate cancer cell lines-all around 48h. Daily percentage of cells proliferating = growth fraction/time between mitoses, where growth fraction is the fraction of cells expressing Ki67 and can be stained for it.
    - Normal prostate cells show a doubling time of about 500 days, with proliferation and death rates evenly matched. This is the steady state under homeostasis.
    - Upon transformation to early stage precancerous neoplastic cells (PINs), proliferation rate increases ~6.9x while death rate inceases ~4x => net growth rate is 0.0045 per day with a doubling time of about 154 days. Late stage PINs show up to 10x higher proliferation and death rates, and therefore no net growth or death. But they show drastically faster turnover (56 days vs 500 days for normal epithelium). Fully transformed cancer cells show a **decrease** of growth and death rates from late stage PINs, but still have net growth compared to normal epithelium at a rate of about 0.0012-0.0014 per day.
    - In both bone and lymph node metastases, there's a big and disproportionate increase in proliferation rates, leading to doubling times of 54 and 33 days respectively.
    - Androgen-independent cancer cells in bone and non-bone sites show similar proliferation rates as the other metastases but higher death rates. They still show net growth at rate of 0.0055-0.0074 per day and doubling times of 125 days and 94 days respectively.
    - See **Table 1: Kinetics parameters of normal and neoplastic prostatic epithelial cells**.

# Lung
1. _Actual growth rate and tumour cell proliferation in human pulmonary neoplasms_ (http://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC1976798&blobtype=pdf)
    - Volumetric doubling times for several histological kinds of lung cancer; shown that cell loss is proportional to LI-high LI also entails high cell loss, but this could also be due to a bad growth model?
    - Well-differentiated tumours tend to have longer DT.
    - **Full list of values in Tables 1 and 2**

# Colon
1. _Accurate reconstruction of the temporal order of mutations in neoplastic progression_ (https://cancerpreventionresearch.aacrjournals.org/content/4/7/1135)
    - Agent-based model of colon crypt used to infer order of kinds of mutations from individual cell-level phlygenies
    - Loss of differentiation is one of the most frequent first mutation event, while second mutation types are more variable.

# Ovary
1. _Differential expression of estrogen receptor subtypes and variants in ovarian cancer_ (http://bmccancer.biomedcentral.com/articles/10.1186/s12885-017-3601-1)
