# =============================================================================
# 单谱图输出_独立版.R - 单谱图输出（全部函数内嵌，单文件可移植，仅 base R）
# 单谱图输出_独立版.R - standalone single-spectrum exporter (self-contained)
#
# 运行：双击/命令行 Rscript 单谱图输出_独立版.R，弹出两个选择窗口：
#   ① 选择 pFind-Filtered.spectra（列：File_Name、Scan_No、Exp.MH+、Charge、
#     Sequence、Modification、Raw_Score、Final_Score、Proteins、Positions）；
#   ② 选择 MGF 文件（可多选）。
# 输出：每谱 16:9 一页 PDF＋300dpi PNG（文件名＝scan<Scan_No>）；
#   内容＝信息行（谱图号/打分/蛋白质|序列|修饰）＋离子梯＋谱图（Base peak/
#   Mass error、0–100 相对强度、同峰合并、8 号子离子原位标注）＋ppm 图。
# 可选环境变量：SELECT=<Scan_No 或行号，逗号分隔>；LIMIT=<最多张数>。
# 依赖：仅 base R；复制到其他电脑可直接运行。
# =============================================================================

# ---- 氨基酸与修饰质量表（内嵌，取自 AA_mass_dataframe.txt）----------------------
# 质量表以文本内嵌于脚本，保证单文件可移植、不依赖外部数据文件
# ---- amino-acid and modification mass table (embedded from AA_mass_dataframe.txt)
# the table is embedded as text so the single file stays portable with no external data（内嵌，取自 AA_mass_dataframe.txt）----------------------
AA_MASSTXT <- "AA	C	H	N	O	S	C12	Mono_weight	N15_weight	C13_weight	
A	3	5	1	1	0	0	71.0371103	72.0341459	74.0471765	
C	3	5	1	1	1	0	103.0091803	104.0062159	106.0192465	
D	4	5	1	3	0	0	115.0269385	116.0239741	119.0403601	
E	5	7	1	3	0	0	129.0425877	130.0396233	134.0593647	
F	9	9	1	1	0	0	147.0684087	148.0654443	156.0986073	
G	2	3	1	1	0	0	57.0214611	58.0184967	59.0281719	
H	6	7	3	1	0	0	137.0589059	140.0500127	143.0790383	
I	6	11	1	1	0	0	113.0840579	114.0810935	119.1041903	
K	6	12	2	1	0	0	128.0949557	130.0890269	134.1150881	
L	6	11	1	1	0	0	113.0840579	114.0810935	119.1041903	
M	5	9	1	1	1	0	131.0404787	132.0375143	136.0572557	
N	4	6	2	2	0	0	114.0429222	116.0369934	118.0563438	
P	5	7	1	1	0	0	97.0527595	98.0497951	102.0695365	
Q	5	8	2	2	0	0	128.0585714	130.0526426	133.0753484	
R	6	12	4	1	0	0	156.1011021	160.0892445	162.1212345	
S	3	5	1	2	0	0	87.0320244	88.02906	90.0420906	
T	4	7	1	2	0	0	101.0476736	102.0447092	105.0610952	
V	5	9	1	1	0	0	99.0684087	100.0654443	104.0851857	
W	11	10	2	1	0	0	186.0793065	188.0733777	197.1162159	
Y	9	9	1	2	0	0	163.0633228	164.0603584	172.0935214
X	6	11	1	1	0	0	113.0840579	114.0810935	119.1041903
U	3	5	1	1	0	0	150.9536356	151.9506705	153.9637001	
Gln->pyro-Glu[AnyN-termQ]	0	-3	-1	0	0	0	-17.026549	-18.023583	-17.026549	H(-3)N(-1)
Glu->pyro-Glu[AnyN-termE]	0	-3	-1	0	0	0	-18.010565	-18.023583	-17.026549	H(-2)O(-1)
2_13C[Any]	2	0	0	0	0	-2	-2.006711	-2.006711	-2.006711	C(2)13C(-2)
1_13C[Any]	1	0	0	0	0	-1	-1.003355	-1.003355	-1.003355	C(1)13C(-1)
Amidated[AnyC-term]	0	1	1	-1	0	0	-0.984016	-0.984016	-0.984016	H(1)N(1)O(-1)
Deamidated[N]	0	-1	-1	1	0	0	0.984016	-0.013019	0.984016	H(-1)N(-1)O(1)
Deamidated[Q]	0	-1	-1	1	0	0	0.984016	-0.013019	0.984016	H(-1)N(-1)O(1)
P_C5N1	5	0	1	0	0	-5	6.013813	6.013813	6.013813	13C(5)C(-5)15N(1)N(-1)
C+12[AnyN-term]	0	0	0	0	0	1	12	12	12	C(1)
Methyl[AnyC-term]	0	2	0	0	0	1	14.01565	14.01565	14.01565	H(2)C(1)
Oxidation[M]	0	0	0	1	0	0	15.994915	15.9949141	15.9949141	O(1)
Cation_Na[AnyC-term]	0	-1	0	0	0	0	21.981943	21.981943	21.981943	H(-1)Na(1)
Cation_Na[D]	0	-1	0	0	0	0	21.981943	21.981943	21.981943	H(-1)Na(1)
Cation_Na[E]	0	-1	0	0	0	0	21.981943	21.981943	21.981943	H(-1)Na(1)
Formyl[AnyN-term]	0	0	0	1	0	1	27.994915	27.994915	27.994915	C(1)O(1)
Acetyl[AnyN-term]	2	2	0	1	0	0	42.010565	42.0105633	44.0172741	H(2)C(2)O(1)
Acetyl[ProteinN-term]	2	2	0	1	0	0	42.010565	42.0105633	44.0172741	H(2)C(2)O(1)
Carbamyl[AnyN-term]	0	1	1	1	0	1	43.005814	43.005814	43.005814	H(1)C(1)N(1)O(1)
Carbamyl[K]	0	1	1	1	0	1	43.005814	43.005814	43.005814	H(1)C(1)N(1)O(1)
Carbamidomethyl[AnyN-term]	0	3	1	1	0	2	57.021464	57.0214611	57.0214611	H(3)C(2)N(1)O(1)
Carbamidomethyl[C]	0	3	1	1	0	2	57.021464	57.0214611	57.0214611	H(3)C(2)N(1)O(1)
Carboxymethyl[C]	0	2	0	2	0	2	58.005479	58.005479	58.005479	H(2)C(2)O(2)
2-hydroxyisobutyrylation[K]	0	6	0	2	0	4	86.036779	86.036779	86.036779	H(6)C(4)O(2)
3-phosphoglyceryl[K]							167.982375	167.982375	167.982375	H(5)C(3)O(6)P(1)
3sulfo[AnyN-term]							183.983029	183.983029	183.983029	H(4)C(7)O(4)S(1)
6[P]	6	0	0	0	0	-6	6	6	6	
AccQTag[AnyN-term]							170.048013	170.048013	170.048013	H(6)C(10)N(2)O(1)
Acetyl[K]							42.010565	42.010565	42.010565	H(2)C(2)O(1)
Acetyl[T]							42.010565	42.010565	42.010565	H(2)C(2)O(1)
Acetyl_13C(2)[ProteinN-term]							44.017274	44.017274	44.017274	H(2)13C(2)O(1)
Acetyl_2H(3)[AnyN-term]							45.029395	45.029395	45.029395	H(-1)2H(3)C(2)O(1)
Acetylhypusine[K]							113.084064	113.084064	113.084064	H(11)C(6)N(1)O(1)
AEBS[H]							183.035399	183.035399	183.035399	H(9)C(8)N(1)O(2)S(1)
AEBS[K]							183.035399	183.035399	183.035399	H(9)C(8)N(1)O(2)S(1)
AEBS[ProteinN-term]							183.035399	183.035399	183.035399	H(9)C(8)N(1)O(2)S(1)
AEBS[Y]							183.035399	183.035399	183.035399	H(9)C(8)N(1)O(2)S(1)
AFB1_Dialdehyde[K]							310.047738	310.047738	310.047738	H(10)C(17)O(6)
Ammonia-loss[N]							-17.026549	-18.023583	-17.026549	H(-3)N(-1)
Ammonium[D]							17.026549	17.026549	17.026549	H(3)N(1)
Arg[AnyN-term]							156.101111	156.101111	156.101111	H(12)C(6)N(4)O(1)
Biotin[AnyN-term]							226.077598	226.077598	226.077598	H(14)C(10)N(2)O(2)S(1)
BITC[AnyN-term]							149.02992	149.02992	149.02992	H(7)C(8)N(1)S(1)
Bromo[W]							77.910511	77.910511	77.910511	H(-1)Br(1)
C8-QAT[AnyN-term]							227.224915	227.224915	227.224915	H(29)C(14)N(1)O(1)
CAF[AnyN-term]							135.983029	135.983029	135.983029	H(4)C(3)O(4)S(1)
Carbamidomethyl[H]							57.021464	57.021464	57.021464	H(3)C(2)N(1)O(1)
Carbamidomethyl[K](Gly[K])							57.021464	57.021464	57.021464	H(3)C(2)N(1)O(1)
Carbamidomethyl[M]							57.021464	57.021464	57.021464	H(3)C(2)N(1)O(1)
Carbamidomethyl[S](Gly[S])							57.021464	57.021464	57.021464	H(3)C(2)N(1)O(1)
Carbofuran[S]							58.029289	58.029289	58.029289	H(4)C(2)N(1)O(1)
Carbonyl[V]							13.979265	13.979265	13.979265	H(-2)O(1)
Carboxy[E]							43.989829	43.989829	43.989829	C(1)O(2)
Carboxyethylpyrrole[K]							122.036779	122.036779	122.036779	H(6)C(7)O(2)
Carboxymethyl[AnyN-term]							58.005479	58.005479	58.005479	H(2)C(2)O(2)
Carboxymethyl[C]							58.005479	58.005479	58.005479	H(2)C(2)O(2)
Cation_Ca[II][D]							37.946941	37.946941	37.946941	H(-2)Ca(1)
Cation_Ca[II][E]							37.946941	37.946941	37.946941	H(-2)Ca(1)
Cation_K[D]							37.955882	37.955882	37.955882	H(-1)K(1)
CHDH[D]							294.183109	294.183109	294.183109	H(26)C(17)O(4)
CLIP_TRAQ_4[AnyN-term]							244.101452	244.101452	244.101452	H(15)C(9)13C(1)N(2)O(5)
Cresylphosphate[T]							170.013281	170.013281	170.013281	H(7)C(7)O(3)P(1)
Cyano[C]							24.995249	24.995249	24.995249	H(-1)C(1)N(1)
Cytopiloyne+water[T]							380.147118	380.147118	380.147118	H(24)C(19)O(8)
Deamidated[R]							0.984016	0.984016	0.984016	H(-1)N(-1)O(1)
Deamidated_18O(1)[N](Delta_H(1)N(-1)18O(1)[N])							2.988261	1.988261	2.988261	H(-1)N(-1)18O(1)
Deamidated_18O(1)[Q]							2.988261	1.988261	2.988261	H(-1)N(-1)18O(1)
Decarboxylation[D]							-30.010565	-30.010565	-30.010565	H(-2)C(-1)O(-1)
Dehydrated[D]							-18.010565	-18.010565	-18.010565	H(-2)O(-1)
Dehydrated[S]							-18.010565	-18.010565	-18.010565	H(-2)O(-1)
Dehydrated[T]							-18.010565	-18.010565	-18.010565	H(-2)O(-1)
Delta_H(2)C(2)[AnyN-term]							26.01565	26.01565	26.01565	H(2)C(2)
Delta_H(3)C(3)O(2)[K]							71.013304	71.013304	71.013304	H(3)C(3)O(2)
Delta_H(-4)O(2)[W]							27.958529	27.958529	27.958529	H(-4)O(2)
Deoxy[T]							-15.994915	-15.994915	-15.994915	O(-1)
DeStreak[C]							75.998285	75.998285	75.998285	H(4)C(2)O(1)S(1)
Dethiomethyl[M]							-48.003371	-48.003371	-48.003371	H(-4)C(-1)S(-1)
dHex[S]							146.057909	146.057909	146.057909	dHex(1)
dHex[T]							146.057909	146.057909	146.057909	dHex(1)
Dicarbamidomethyl[AnyN-term]							114.042927	114.042927	114.042927	H(6)C(4)N(2)O(2)
Didehydro[S]							-2.01565	-2.01565	-2.01565	H(-2)
Diethyl[AnyN-term]							56.0626	56.0626	56.0626	H(8)C(4)
Diethylphosphothione[H]							152.006087	152.006087	152.006087	H(9)C(4)O(2)P(1)S(1)
Diisopropylphosphate[K]							164.060231	164.060231	164.060231	H(13)C(6)O(3)P(1)
Dimethyl[AnyN-term](Ethyl[AnyN-term])							28.0313	28.0313	28.0313	H(4)C(2)
Dimethyl[N]							28.0313	28.0313	28.0313	H(4)C(2)
Dimethyl_2H(4)[AnyN-term]							32.056407	32.056407	32.056407	2H(4)C(2)
Dimethyl_2H(6)[AnyN-term]							34.068961	34.068961	34.068961	H(-2)2H(6)C(2)
Dioxidation[E]							31.989829	31.989829	31.989829	O(2)
Dioxidation[I]							31.989829	31.989829	31.989829	O(2)
Dioxidation[M]							31.989829	31.989829	31.989829	O(2)
Dioxidation[W]							31.989829	31.989829	31.989829	O(2)
EDT-iodoacetyl-PEO-biotin[S]							490.174218	490.174218	490.174218	H(34)C(20)N(4)O(4)S(3)
EDT-iodoacetyl-PEO-biotin[T]							490.174218	490.174218	490.174218	H(34)C(20)N(4)O(4)S(3)
Ethanolyl[K]							44.026215	44.026215	44.026215	H(4)C(2)O(1)
Ethylphosphate[K]							107.997631	107.997631	107.997631	H(5)C(2)O(3)P(1)
Galactosyl[AnyN-term]							178.047738	178.047738	178.047738	O(1)Hex(1)
GG[C](Dicarbamidomethyl[C])							114.042927	114.042927	114.042927	H(6)C(4)N(2)O(2)
GG[K](Dicarbamidomethyl[K])							114.042927	114.042927	114.042927	H(6)C(4)N(2)O(2)
GGQ[K]							242.101505	242.101505	242.101505	H(14)C(9)N(4)O(4)
GIST-Quat[AnyN-term]							127.099714	127.099714	127.099714	H(13)C(7)N(1)O(1)
GIST-Quat[K]							127.099714	127.099714	127.099714	H(13)C(7)N(1)O(1)
GIST-Quat_2H(3)[AnyN-term]							130.118544	130.118544	130.118544	H(10)2H(3)C(7)N(1)O(1)
GIST-Quat_2H(3)[K]							130.118544	130.118544	130.118544	H(10)2H(3)C(7)N(1)O(1)
glycidamide[AnyN-term]							87.032028	87.032028	87.032028	H(5)C(3)N(1)O(2)
Guanidinyl[AnyN-term]							42.021798	42.021798	42.021798	H(2)C(1)N(2)
Hex(1)HexA(1)[T]							338.084912	338.084912	338.084912	Hex(1)HexA(1)
Hex(1)Pent(2)[S]							426.137341	426.137341	426.137341	Pent(2)Hex(1)
Hex[AnyN-term]							162.052824	162.052824	162.052824	Hex(1)
Hex[T]							162.052824	162.052824	162.052824	Hex(1)
HexN[T]							161.068808	161.068808	161.068808	HexN(1)
HexNAc(1)dHex(1)[T]							349.137281	349.137281	349.137281	dHex(1)HexNAc(1)
HexNAc(1)NeuAc(1)[S]							494.174789	494.174789	494.174789	HexNAc(1)NeuAc(1)
HexNAc[T]							203.079373	203.079373	203.079373	HexNAc(1)
HNE[A]							156.11503	156.11503	156.11503	H(16)C(9)O(2)
HNE[L]							156.11503	156.11503	156.11503	H(16)C(9)O(2)
HNE-Delta_H(2)O[H]							138.104465	138.104465	138.104465	H(14)C(9)O(1)
Homocysteic_acid[M]							33.969094	33.969094	33.969094	H(-2)C(-1)O(3)
ICPL_13C(6)[AnyN-term]							111.041593	111.041593	111.041593	H(3)13C(6)N(1)O(1)
ICPL_13C(6)2H(4)[AnyN-term]							115.0667	115.0667	115.0667	H(-1)2H(4)13C(6)N(1)O(1)
ICPL_2H(4)[AnyN-term]							109.046571	109.046571	109.046571	H(-1)2H(4)C(6)N(1)O(1)
Iminobiotin[AnyN-term]							225.093583	225.093583	225.093583	H(15)C(10)N(3)O(1)S(1)
Iodoacetanilide[AnyN-term]							133.052764	133.052764	133.052764	H(7)C(8)N(1)O(1)
Iodoacetanilide_13C(6)[AnyN-term]							139.072893	139.072893	139.072893	H(7)C(2)13C(6)N(1)O(1)
iodoTMT[H]							324.216141	324.216141	324.216141	H(28)C(16)N(4)O(3)
ISD_z+2_ion[AnyN-term]							-15.010899	-16.010899	-15.010899	H(-1)N(-1)
iTRAQ4plex114[AnyN-term]							144.105918	144.105918	144.105918	H(12)C(5)13C(2)N(2)18O(1)
Lys[AnyN-term]							128.094963	128.094963	128.094963	H(12)C(6)N(2)O(1)
Malonyl[S]							86.000394	86.000394	86.000394	H(2)C(3)O(3)
Menadione-HQ[C]							172.05243	172.05243	172.05243	H(8)C(11)O(2)
MercaptoEthanol[S]							60.003371	60.003371	60.003371	H(4)C(2)S(1)
Methyl_2H(3)13C(1)[AnyN-term]							18.037835	18.037835	18.037835	H(-1)2H(3)13C(1)
Methylphosphonate[T]							77.987066	77.987066	77.987066	H(3)C(1)O(2)P(1)
Methylthio[AnyN-term]							45.987721	45.987721	45.987721	H(2)C(1)S(1)
Met-loss+Acetyl[ProteinN-termM]							-89.02992	-90.02992	-89.02992	H(-7)C(-3)N(-1)S(-1)
mTRAQ[AnyN-term]							140.094963	140.094963	140.094963	H(12)C(7)N(2)O(1)
mTRAQ_13C(6)15N(2)[AnyN-term]							148.109162	148.109162	148.109162	H(12)C(1)13C(6)15N(2)O(1)
Myristoyl[C]							210.198366	210.198366	210.198366	H(26)C(14)O(1)
NDA[AnyN-term]							175.042199	175.042199	175.042199	H(5)C(13)N(1)
NDA[K]							175.042199	175.042199	175.042199	H(5)C(13)N(1)
NHS-LC-Biotin[AnyN-term]							339.161662	339.161662	339.161662	H(25)C(16)N(3)O(3)S(1)
NIC[AnyN-term](ICPL[AnyN-term])							105.021464	105.021464	105.021464	H(3)C(6)N(1)O(1)
Octanoyl[S]							126.104465	126.104465	126.104465	H(14)C(8)O(1)
O-Methylphosphate[T]							93.981981	93.981981	93.981981	H(3)C(1)O(3)P(1)
O-pinacolylmethylphosphonate[S]							162.080967	162.080967	162.080967	H(15)C(7)O(2)P(1)
Oxidation[S]							15.994915	15.994915	15.994915	O(1)
Oxidation[W]							15.994915	15.994915	15.994915	O(1)
Pent(2)[S]							264.084518	264.084518	264.084518	Pent(2)
PET[T]							121.035005	121.035005	121.035005	H(7)C(7)N(1)O(-1)S(1)
Phenylisocyanate_2H(5)[AnyN-term]							124.068498	124.068498	124.068498	2H(5)C(7)N(1)O(1)
Phospho[E]							79.966331	79.966331	79.966331	H(1)O(3)P(1)
Phosphoadenosine[K]							329.05252	329.05252	329.05252	H(12)C(10)N(5)O(6)P(1)
Phosphoadenosine[T]							329.05252	329.05252	329.05252	H(12)C(10)N(5)O(6)P(1)
Phosphopantetheine[S]							340.085794	340.085794	340.085794	H(21)C(11)N(2)O(6)P(1)S(1)
Propionamide[AnyN-term]							71.037114	71.037114	71.037114	H(5)C(3)N(1)O(1)
Propionamide[C]							71.037114	71.037114	71.037114	H(5)C(3)N(1)O(1)
Propionyl[AnyN-term]							56.026215	56.026215	56.026215	H(4)C(3)O(1)
Propionyl_13C(3)[AnyN-term]							59.036279	59.036279	59.036279	H(4)13C(3)O(1)
Pyridylacetyl[AnyN-term](Phenylisocyanate[AnyN-term])							119.037114	119.037114	119.037114	H(5)C(7)N(1)O(1)
pyrophospho[T]							159.932662	159.932662	159.932662	H(2)O(6)P(2)
Quinone[W]							29.974179	29.974179	29.974179	H(-2)O(2)
R004							3.988142	3.988142	3.988142	N(-4)15N(4)
Sulfide[D]							31.972071	31.972071	31.972071	S(1)
sulfo+amino[Y]							94.967714	94.967714	94.967714	H(1)N(1)O(3)S(1)
Thiazolidine[W]							12	12	12	C(1)
TMAB[AnyN-term]							128.107539	128.107539	128.107539	H(14)C(7)N(1)O(1)
TMT2plex[AnyN-term]							225.155833	225.155833	225.155833	H(20)C(11)13C(1)N(2)O(2)
trifluoro[L]							53.971735	53.971735	53.971735	H(-3)F(3)
Trioxidation[C]							47.984744	47.984744	47.984744	O(3)
Trioxidation[W]							47.984744	47.984744	47.984744	O(3)
UgiJoullieProGly[D]							154.074228	154.074228	154.074228	H(10)C(7)N(2)O(2)
Unknown_248[AnyN-term]							248.19876	248.19876	248.19876	H(28)C(13)O(4)
Unknown_306[AnyN-term]							306.095082	306.095082	306.095082	H(18)C(12)O(9)
Unknown_306[D]							306.095082	306.095082	306.095082	H(18)C(12)O(9)
Unknown_420[AnyN-term]							420.050846	420.050846	420.050846	H(14)C(10)N(9)O(6)S(2)
Xlink_BS2G[114][K](Gluratylation[K])							114.031694	114.031694	114.031694	H(6)C(5)O(3)
Xlink_BS2G[96][K]							96.021129	96.021129	96.021129	H(4)C(5)O(2)
Xlink_BuUrBu[196][K]							196.084792	196.084792	196.084792	H(12)C(9)N(2)O(3)
Xlink_DSS[156][ProteinN-term]							156.078644	156.078644	156.078644	H(12)C(8)O(3)
Xlink_DSSO[54][K](Delta_H(2)C(3)O(1)[K])							54.010565	54.010565	54.010565	H(2)C(3)O(1)
Xlink_DST[132][K]							132.005873	132.005873	132.005873	H(4)C(4)O(5)
Xlink_DST[56][ProteinN-term]							55.989829	55.989829	55.989829	C(2)O(2)
Xlink_DTSSP[88][K]							87.998285	87.998285	87.998285	H(4)C(3)O(1)S(1)
ZQG[K]							320.100836	320.100836	320.100836	H(16)C(15)N(2)O(6)、
Label_13C(6)[K]							6.020129	6.020129	6.020129	C(-6)13C(6)
Label_13C(6)15N(4)[R]							10.008269	10.008269	10.008269	C(-6)13C(6)N(-4)15N(4)"

# ---- 基础小工具（base R 实现的原 compat.R 同名函数）----------------------------
# 三个字符串小函数以 base R 重实现，替代 stringr 依赖
# ---- small base-R helpers (the same names as the original compat.R layer)
# three string helpers reimplemented in base R to replace the stringr dependency--------------------------------------------------------
str_detect <- function(string, pattern) grepl(pattern, string, fixed = TRUE)
str_split  <- function(string, pattern) strsplit(string, pattern, fixed = TRUE)

# ---- MGF 读取（向量化解析＋并行分块，仅 base R parallel）------------------------
# 单文件内实现：S4 类定义内嵌、逐文件 readLines 后按 BEGIN/END 切块，
# 谱数多时以 PSOCK 集群分块并行解析（clusterExport 直接传函数，worker 免 source）
# ---- MGF reading (vectorized parse + parallel chunks, base parallel only)
# self-contained: the S4 class is defined inline, each file is read once and
# split on BEGIN/END; large files are chunk-parsed on a PSOCK cluster (the
# worker receives the parser via clusterExport, so no sourcing is needed)（向量化解析＋并行分块，base R parallel）----------------------------
MGF_CLASS_DEF <- paste0(
  "setClass('MGF', ",
  "slots = list(title = 'character', charge = 'integer', ",
  "rtInSeconds = 'numeric', pepMass = 'numeric', spectrum = 'data.frame'), ",
  "prototype = list(title = character(0), charge = integer(0), ",
  "rtInSeconds = numeric(0), pepMass = numeric(0), ",
  "spectrum = data.frame()))")
eval(parse(text = MGF_CLASS_DEF))

# 块级向量化解析（原库 parse_lines_vectorized 逐字移植）：一次正则分类、
# 两次 cumsum、一次 strsplit＋as.numeric 解析全部峰行；返回 S4 MGF 对象列表。
parse_lines_vectorized <- function(lines) {
  if (length(lines) == 0) return(list())
  is_peak <- grepl("^[0-9.]", lines)
  is_beg  <- startsWith(lines, "BEGIN IONS")
  spec_id <- cumsum(is_beg)
  in_spec <- spec_id > 0L
  max_spectra <- sum(is_beg)
  pk    <- lines[is_peak & in_spec]
  pk_id <- spec_id[is_peak & in_spec]
  n_peaks <- length(pk)
  cnt <- tabulate(pk_id, nbins = max_spectra)
  peak_ends <- as.integer(cumsum(cnt))
  hdr_val <- function(prefix) {
    i <- which(in_spec & startsWith(lines, prefix))
    i <- i[!duplicated(spec_id[i], fromLast = TRUE)]
    v <- rep(NA_character_, max_spectra)
    v[spec_id[i]] <- sub(prefix, "", lines[i], fixed = TRUE)
    v
  }
  titles  <- hdr_val("TITLE="); titles[is.na(titles)] <- ""
  ch_raw  <- hdr_val("CHARGE=")
  charges <- rep(NA_integer_, max_spectra)
  ch_ok   <- grepl("^\\d+\\+$", ch_raw)
  charges[ch_ok] <- as.integer(sub("\\+$", "", ch_raw[ch_ok]))
  rt_raw  <- hdr_val("RTINSECONDS=")
  rts     <- ifelse(is.na(rt_raw), 0, suppressWarnings(as.numeric(rt_raw)))
  pms     <- suppressWarnings(as.numeric(sub("[ \t]+.*$", "", hdr_val("PEPMASS="))))
  if (n_peaks > 0L) {
    pairs <- strsplit(pk, "[ \t]+")
    nums <- as.numeric(unlist(pairs, use.names = FALSE))
    if (length(nums) == 2L * n_peaks) {
      mz_all  <- nums[seq.int(1L, length(nums), 2L)]
      int_all <- nums[seq.int(2L, length(nums), 2L)]
      rm(pairs, nums)
    } else {
      mz_all  <- vapply(pairs, function(p) as.numeric(p[1]), 0.0)
      int_all <- vapply(pairs, function(p) as.numeric(p[2]), 0.0)
    }
  } else {
    mz_all <- numeric(0)
    int_all <- numeric(0)
  }
  keep <- which(cnt > 0L)
  lo <- c(0L, peak_ends[-length(peak_ends)])[keep] + 1L
  hi <- peak_ends[keep]
  out <- vector("list", length(keep))
  for (q in seq_along(keep)) {
    seg <- seq.int(lo[q], hi[q])
    out[[q]] <- new("MGF",
                    title = titles[keep[q]],
                    charge = charges[keep[q]],
                    rtInSeconds = rts[keep[q]],
                    pepMass = pms[keep[q]],
                    spectrum = structure(
                      list(m_z = mz_all[seg], Intensity = int_all[seg]),
                      class = "data.frame", row.names = .set_row_names(length(seg))))
  }
  out
}

# 逐文件读取：谱数多时按块并行解析（clusterExport 内嵌函数，worker 无需 source）
read_mgf <- function(files) {
  titles_all <- list(); objs_all <- list()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    i_beg <- which(lines == "BEGIN IONS"); i_end <- which(lines == "END IONS")
    if (!length(i_beg)) next
    raw_name <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(f)))
    n_spec <- length(i_beg)
    ncores <- parallel::detectCores()
    if (is.na(ncores) || ncores < 2L) ncores <- 2L
    ncores <- max(1L, min(8L, ncores - 1L))
    chunk_size <- 800L
    chunks <- split(seq_len(n_spec), ceiling(seq_len(n_spec) / chunk_size))
    if (length(chunks) > 1L && ncores > 1L) {
      cl <- parallel::makeCluster(min(ncores, length(chunks)))
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("parse_lines_vectorized", "MGF_CLASS_DEF"),
                              envir = environment())
      parallel::clusterEvalQ(cl, eval(parse(text = MGF_CLASS_DEF)))
      res <- parallel::parLapply(cl, chunks, function(ci) {
        lo <- i_beg[ci[1]]; hi <- i_end[tail(ci, 1)]
        parse_lines_vectorized(lines[lo:hi])
      })
      objs <- do.call(c, res)
    } else {
      objs <- parse_lines_vectorized(lines[i_beg[1]:i_end[length(i_end)]])
    }
    tt <- vapply(objs, function(o) o@title, "")
    oo <- lapply(objs, function(o) list(title = o@title, charge = o@charge,
                                        pepMass = o@pepMass, spectrum = o@spectrum))
    titles_all[[raw_name]] <- tt; objs_all[[raw_name]] <- oo
    cat(sprintf("读取 %s：%d 谱（%d 块/%d 核）
", basename(f), n_spec,
                length(chunks), ncores))
  }
  list(objs = objs_all, titles = titles_all)
}

title_raw_scan <- function(title) {
  p <- strsplit(title, ".", fixed = TRUE)[[1]]
  c(raw = p[1], scan = p[2])
}
title_scan_key <- function(title) {
  p <- strsplit(title, ".", fixed = TRUE)[[1]]
  paste0(p[1], ".", p[2])
}
find_spec <- function(lu, fn) {
  raw <- unname(title_raw_scan(fn)["raw"])
  objs <- lu$objs[[raw]]
  if (is.null(objs)) return(NULL)
  li <- match(fn, lu$titles[[raw]])
  if (is.na(li)) {
    key <- title_scan_key(fn)
    ci <- which(startsWith(lu$titles[[raw]], paste0(key, ".")))
    if (length(ci) == 0) return(NULL)
    li <- ci[1]
  }
  objs[[li]]
}

# ---- 氨基酸对照表（三字母/单字母/全名，原库同款）-------------------------------
# 仅用于 X->Y 替换修饰的三字母名到单字母码转换
# ---- amino-acid lookup table (three-letter / one-letter / full name)
# used only to convert the three-letter names of X->Y substitution mods（原库同款）----------------------------------------------------
amino_acids_three_letter <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
                              "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp",
                              "Tyr","Val","Xle")
amino_acids_one_letter <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F",
                            "P","S","T","W","Y","V","I")
amino_acid_names <- c("Alanine","Arginine","Asparagine","Aspartic acid","Cysteine",
                      "Glutamine","Glutamic acid","Glycine","Histidine","Isoleucine",
                      "Leucine","Lysine","Methionine","Phenylalanine","Proline",
                      "Serine","Threonine","Tryptophan","Tyrosine","Valine",
                      "Isoleucine_or_Leucine")
amino_acid_table <- data.frame(Full_Name = amino_acid_names,
                               Three_Letter = amino_acids_three_letter,
                               One_Letter = amino_acids_one_letter)
rm(amino_acids_three_letter, amino_acids_one_letter, amino_acid_names)

# ---- 理论 b/y 离子计算（原 by_ions_position 的 Unlabel 分支逐字移植）-------------
# 输出 ions_table（b/y 1+ 与 2+ 质量）与 pep_mass（含质子的单同位素质量）
# ---- theoretical b/y ion computation (verbatim port of the Unlabel branch)
# returns the ions_table (b/y 1+ and 2+ masses) and pep_mass (monoisotopic
# mass including the proton)（原 by_ions_position 的 Unlabel 分支逐字移植）------------------
by_ions_position <- function(peptide = "", modification = "", label = "Unlabel",
                             AA_mass_dataframe, amino_acid_table) {
  if (nchar(peptide) == 0) return("Not a peptide.")
  H_atom_weight <- 1.007825
  O_atom_weight <- 15.994914
  N_terminal_weight <- H_atom_weight
  C_terminal_weight <- H_atom_weight + O_atom_weight
  AA_weight <- AA_mass_dataframe$Mono_weight[
    match(strsplit(peptide, "", fixed = TRUE)[[1]], AA_mass_dataframe$AA)]
  if (modification != "") {
    modification_split <- str_split(modification, ";")[[1]]
    modification_split <- modification_split[modification_split != ""]
    for (modification_index in seq_along(modification_split)) {
      temp_modification_split <- str_split(modification_split[modification_index], ",")[[1]]
      if (!temp_modification_split[2] %in% AA_mass_dataframe$AA) {
        if (str_detect(temp_modification_split[2], "->")) {
          AA_to_AA <- str_split(temp_modification_split[2], "->")[[1]]
          AA_changed <- strsplit(AA_to_AA[2], "[", fixed = TRUE)[[1]][1]
          AA_changed <- amino_acid_table$One_Letter[
            which(amino_acid_table$Three_Letter == AA_changed)]
          if (length(AA_changed) != 1) return(NA)
          AA_weight[as.numeric(temp_modification_split[1])] <-
            AA_mass_dataframe$Mono_weight[
              which(AA_mass_dataframe$AA == AA_changed)]
        } else return(NA)
      } else {
        mass_shift <- AA_mass_dataframe$Mono_weight[
          which(AA_mass_dataframe$AA == temp_modification_split[2])]
        if (str_detect(temp_modification_split[2], "N-term")) {
          N_terminal_weight <- N_terminal_weight + mass_shift
        } else if (str_detect(temp_modification_split[2], "C-term")) {
          C_terminal_weight <- C_terminal_weight + mass_shift
        } else {
          AA_weight[as.numeric(temp_modification_split[1])] <-
            AA_weight[as.numeric(temp_modification_split[1])] + mass_shift
        }
      }
    }
  }
  ions_table <- as.data.frame(matrix(NA, nrow = max(0L, nchar(peptide) - 1L), ncol = 5))
  colnames(ions_table) <- c("Number", "b_1_charge", "b_2_charge",
                            "y_1_charge", "y_2_charge")
  ions_table$Number <- seq_len(max(0L, nchar(peptide) - 1L))
  ions_table$b_1_charge <- cumsum(AA_weight)[seq_len(nchar(peptide) - 1L)] +
    N_terminal_weight
  ions_table$y_1_charge <- cumsum(rev(AA_weight))[seq_len(nchar(peptide) - 1L)] +
    C_terminal_weight
  ions_table$y_1_charge <- ions_table$y_1_charge + H_atom_weight * 2
  ions_table$b_2_charge <- ions_table$b_1_charge / 2 + H_atom_weight / 2
  ions_table$y_2_charge <- ions_table$y_1_charge / 2 + H_atom_weight / 2
  pep_mass <- sum(AA_weight, N_terminal_weight, C_terminal_weight) + H_atom_weight
  list(ions_table = ions_table, pep_mass = pep_mass)
}

# ---- 窗口强度与覆盖计数（原库同款，排序＋前缀和，O(n log n)）-------------------
# win_intensity＝窗口内强度之和；cover_count＝每个点被多少窗口覆盖
# ---- window intensity and coverage counts (library originals; sort + prefix
# sums, O(n log n)) - win_intensity sums the intensity inside each window and
# cover_count counts the windows covering each point（原库同款）------------------------------------------------
win_intensity <- function(lo, hi, mz, int) {
  o <- order(mz)
  mz_s <- mz[o]
  ps <- c(0, cumsum(int[o]))
  n_lt <- findInterval(lo - 1e-9, mz_s)
  n_le <- findInterval(hi, mz_s)
  out <- ps[n_le + 1L] - ps[n_lt + 1L]
  out[is.na(lo) | is.na(hi)] <- 0
  out
}
cover_count <- function(px, lo, hi) {
  ok <- !is.na(lo) & !is.na(hi)
  lo <- lo[ok]; hi <- hi[ok]
  findInterval(px, sort(lo)) - findInterval(px - 1e-9, sort(hi))
}

# ---- b/y 离子匹配与打分（原 ions_match 逐字移植；melt 以直接长格式代替）---------
# PSM_score_by＝b/y 匹配离子强度占总强度百分比；PSM_score_all＝含 a/水/氨损失
# ---- b/y ion matching and scoring (verbatim port of ions_match; the melt step
# is replaced by a direct long-format build) - PSM_score_by is the share of the
# total intensity carried by matched b/y ions; PSM_score_all also counts a ions
# and the water/ammonia losses（原 ions_match 逐字移植；melt 以直接长格式代替）-----------
ions_match <- function(therotical_by_ions, mgf_ions,
                       mass_shift_tolerence = 20, ppm_or_Da = c("ppm", "Da")) {
  ions_table <- therotical_by_ions$ions_table
  by_ions_long <- data.frame(
    Type = rep(colnames(ions_table)[2:5], each = nrow(ions_table)),
    m_z = as.numeric(unlist(ions_table[, 2:5])),
    stringsAsFactors = FALSE)
  by_ions_long$Charge <- 0
  by_ions_long$Charge[by_ions_long$Type == "b_1_charge"] <- 1
  by_ions_long$Charge[by_ions_long$Type == "y_1_charge"] <- 1
  by_ions_long$Charge[by_ions_long$Type == "b_2_charge"] <- 2
  by_ions_long$Charge[by_ions_long$Type == "y_2_charge"] <- 2
  by_ions_long$Type[by_ions_long$Type %in% c("b_1_charge", "b_2_charge")] <- "b"
  by_ions_long$Type[by_ions_long$Type %in% c("y_1_charge", "y_2_charge")] <- "y"
  H2O_loss <- 18.010565
  NH3_loss <- 17.026549
  CO_loss <- 27.994915
  ions_H2O_loss <- by_ions_long
  ions_H2O_loss$m_z <- ions_H2O_loss$m_z - H2O_loss / ions_H2O_loss$Charge
  ions_NH3_loss <- by_ions_long
  ions_NH3_loss$m_z <- ions_NH3_loss$m_z - NH3_loss / ions_NH3_loss$Charge
  a_ions <- by_ions_long[by_ions_long$Type == "b", ]
  a_ions$m_z <- a_ions$m_z - CO_loss / a_ions$Charge
  combined_ions <- rbind(by_ions_long, ions_H2O_loss, ions_NH3_loss, a_ions)
  ppm_or_Da <- tolower(ppm_or_Da[1])
  if (!ppm_or_Da %in% c("ppm", "da") || !is.numeric(mass_shift_tolerence))
    return(NA)
  if (ppm_or_Da == "ppm") {
    by_ions_long$low_m_z <- by_ions_long$m_z * (1 - mass_shift_tolerence * 1e-6)
    by_ions_long$high_m_z <- by_ions_long$m_z * (1 + mass_shift_tolerence * 1e-6)
    combined_ions$low_m_z <- combined_ions$m_z * (1 - mass_shift_tolerence * 1e-6)
    combined_ions$high_m_z <- combined_ions$m_z * (1 + mass_shift_tolerence * 1e-6)
  } else {
    by_ions_long$low_m_z <- by_ions_long$m_z - mass_shift_tolerence / by_ions_long$Charge
    by_ions_long$high_m_z <- by_ions_long$m_z + mass_shift_tolerence / by_ions_long$Charge
    combined_ions$low_m_z <- combined_ions$m_z - mass_shift_tolerence / combined_ions$Charge
    combined_ions$high_m_z <- combined_ions$m_z + mass_shift_tolerence / combined_ions$Charge
  }
  by_ions_long$match_intensity <- win_intensity(by_ions_long$low_m_z,
    by_ions_long$high_m_z, mgf_ions$m_z, mgf_ions$Intensity)
  mgf_ions$Match <- cover_count(mgf_ions$m_z, by_ions_long$low_m_z,
                                by_ions_long$high_m_z)
  mgf_ions$All_Match <- cover_count(mgf_ions$m_z, combined_ions$low_m_z,
                                    combined_ions$high_m_z)
  by_ions_coverage <- round(length(which(by_ions_long$match_intensity > 0)) /
                            nrow(by_ions_long) * 100, 2)
  PSM_score_by <- round(sum(mgf_ions$Intensity[which(mgf_ions$Match > 0)]) /
                        sum(mgf_ions$Intensity) * 100, 2)
  PSM_score_all <- round(sum(mgf_ions$Intensity[which(mgf_ions$All_Match > 0)]) /
                         sum(mgf_ions$Intensity) * 100, 2)
  c(by_ions_long = by_ions_long, by_ions_coverage = by_ions_coverage,
    PSM_score_by = PSM_score_by, PSM_score_all = PSM_score_all)
}

# ---- 输入：交互式文件选择（choose.files 弹窗；不可用时退回手动输入）------------
# ① 选择 pFind-Filtered.spectra；② 选择 MGF 文件（可多选，文件夹自动展开 .mgf）
# ---- input: interactive file selection (choose.files dialogs; manual entry as
# the fallback) - 1) choose pFind-Filtered.spectra; 2) choose the MGF files
# (multi-select; folders are expanded to their .mgf files) --------------------------------------------------------
cat("请选择 pFind-Filtered.spectra 文件...\n")
PF_FILE <- choose.files(default = "", caption = "选择 pFind-Filtered.spectra",
                        multi = FALSE)
if (!length(PF_FILE) || !nzchar(PF_FILE)) stop("未选择 pFind-Filtered.spectra")
cat("请选择 MGF 文件（可多选）...\n")
MGF_FILES <- choose.files(default = "", caption = "选择 MGF 文件（可多选）",
                          multi = TRUE)
if (!length(MGF_FILES)) {
  cat("未选择 MGF：可输入 MGF 文件或文件夹路径（逗号分隔）：")
  inp <- readline()
  if (!nzchar(inp)) stop("未提供 MGF")
  MGF_FILES <- trimws(strsplit(inp, "[,;]")[[1]])
}
for (p in MGF_FILES) {
  if (dir.exists(p))
    MGF_FILES <- c(MGF_FILES, list.files(p, "\\.mgf$", full.names = TRUE))
}
MGF_FILES <- unique(MGF_FILES[file.exists(MGF_FILES) & grepl("\\.mgf$", MGF_FILES)])

pf <- read.delim(PF_FILE, sep = "\t", check.names = FALSE,
                 na.strings = character(0), quote = "", fill = TRUE)
need_cols <- c("File_Name", "Scan_No", "Exp.MH+", "Charge", "Sequence", "Modification")
miss_c <- setdiff(need_cols, names(pf))
if (length(miss_c)) stop("pFind-Filtered.spectra 缺列：", paste(miss_c, collapse = ","))
for (cn in need_cols) pf[[cn]] <- as.character(pf[[cn]])
pf$Scan_No <- suppressWarnings(as.numeric(pf$Scan_No))
pf$Charge <- suppressWarnings(as.integer(pf$Charge))
pf$ExpMH <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pf[["Exp.MH+"]])))
if (!"Raw_Score" %in% names(pf)) pf$Raw_Score <- "-"
if (!"Final_Score" %in% names(pf)) pf$Final_Score <- "-"
if (!"Proteins" %in% names(pf)) pf$Proteins <- ""
if (!"Positions" %in% names(pf)) pf$Positions <- ""

SEL_ENV <- Sys.getenv("SELECT", "")
LIMIT <- as.numeric(Sys.getenv("LIMIT", "Inf"))
if (nzchar(SEL_ENV)) {
  want <- suppressWarnings(as.numeric(trimws(strsplit(SEL_ENV, "[,;]")[[1]])))
  sel <- which(pf$Scan_No %in% want | seq_len(nrow(pf)) %in% want)
} else sel <- seq_len(nrow(pf))
if (is.finite(LIMIT)) sel <- head(sel, LIMIT)
cat(sprintf("pFind-Filtered %d 行；选择 %d 行\n", nrow(pf), length(sel)))

aa <- read.delim(text = AA_MASSTXT, sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table
# 多 MGF：pFind-Filtered 按 MGF（raw）拆分，逐 MGF 读取与输出（子文件夹＝MGF 文件茎）
mgf_stems <- sub("_HCDFT$", "", sub("\\.mgf$", "", basename(MGF_FILES)))
pf_raw <- vapply(strsplit(as.character(pf$File_Name), ".", fixed = TRUE),
                 function(p) p[1], "")
pf_mgf_idx <- match(pf_raw, mgf_stems)                    # 行 → MGF 文件序号

# ---- 绘制引擎（与 Fig.7 谱图定稿版式一致）--------------------------------------
# 颜色表与优先序、理论离子梯（a/b/c/x/y/z 1+/2+ 与 M 系）、20 ppm 窗匹配、
# 同峰合并标签（b 优先）、离子梯（黑色 1 磅切口竖条、价态前缀、修饰橙色加粗）
# ---- the drawing engine (identical to the finalized Fig.7 spectrum layout)
# the color table and priority order, the theoretical ion ladder (a/b/c/x/y/z
# 1+/2+ plus the M series), 20 ppm window matching, merged same-peak labels
# (b first), and the ion ladder (black 1-pt ticks, charge prefix, orange bold
# modified residues)（与谱图可视化同格式）------------------------------------------------
ION_COL <- c(a = "#004000", b = "#008000", c = "#CA1D52", x = "#650E29",
             y = "#CA1D52", z = "#943AA4", M = "#808080", I = "#FF8C00", O = "#9ACD32")
ION_PRI <- c(b = 1, y = 2, a = 3, c = 4, x = 5, z = 6, M = 7, I = 8, O = 9)
CO_M <- 27.994915; NH3_M <- 17.026549; H2O_M <- 18.010565; H_M <- 1.007276

all_ions <- function(seq, mod, mz, ints, tol_ppm = 20) {
  for (md in unique(c(mod, ""))) {
    bip <- NULL
    tt <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(seq, md, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(tt, "try-error") || is.null(bip) || length(bip) == 1 ||
        is.null(bip$ions_table)) next
    t <- bip$ions_table; n <- t$Number
    theo <- data.frame(
      mz  = c(t$b_1_charge, t$b_1_charge - CO_M, t$b_1_charge + NH3_M,
              t$b_2_charge, t$b_2_charge - CO_M/2, t$b_2_charge + NH3_M/2,
              t$y_1_charge, t$y_1_charge + CO_M, t$y_1_charge - NH3_M,
              t$y_2_charge, t$y_2_charge + CO_M/2, t$y_2_charge - NH3_M/2),
      lab = c(paste0("b", n, "+"), paste0("a", n, "+"), paste0("c", n, "+"),
              paste0("b", n, "++"), paste0("a", n, "++"), paste0("c", n, "++"),
              paste0("y", n, "+"), paste0("x", n, "+"), paste0("z", n, "+"),
              paste0("y", n, "++"), paste0("x", n, "++"), paste0("z", n, "++")),
      ion = c(rep(c("b", "a", "c"), each = length(n)),
              rep(c("b", "a", "c"), each = length(n)),
              rep(c("y", "x", "z"), each = length(n)),
              rep(c("y", "x", "z"), each = length(n))),
      num = c(rep(n, 3), rep(n, 3), rep(n, 3), rep(n, 3)),
      chg = c(rep(1L, 3 * length(n)), rep(2L, 3 * length(n)),
              rep(1L, 3 * length(n)), rep(2L, 3 * length(n))),
      stringsAsFactors = FALSE)
    if (is.finite(bip$pep_mass)) {
      pm <- bip$pep_mass
      theo <- rbind(theo, data.frame(
        mz = c(pm + H_M, pm + H_M - H2O_M, pm + H_M - NH3_M,
               (pm + 2 * H_M) / 2, (pm + 2 * H_M - H2O_M) / 2, (pm + 2 * H_M - NH3_M) / 2),
        lab = c("M+", "M-H2O+", "M-NH3+", "M++", "M-H2O++", "M-NH3++"),
        ion = "M", num = NA_integer_, chg = c(1L, 1L, 1L, 2L, 2L, 2L),
        stringsAsFactors = FALSE))
    }
    theo <- theo[is.finite(theo$mz), ]
    lo <- theo$mz * (1 - tol_ppm * 1e-6); hi <- theo$mz * (1 + tol_ppm * 1e-6)
    i_lo <- findInterval(lo, mz); i_hi <- findInterval(hi, mz)
    peak_i <- rep(NA_integer_, nrow(theo))
    for (j in seq_len(nrow(theo))) {
      if (i_hi[j] <= i_lo[j]) next
      cand <- (i_lo[j] + 1L):i_hi[j]
      if (length(cand)) peak_i[j] <- cand[which.max(ints[cand])]
    }
    theo$peak <- peak_i
    return(theo)
  }
  NULL
}

draw_ladder <- function(sq, ions, mod = "", chg_v = 0L) {
  n <- nchar(sq)
  if (n < 2L) return(invisible())
  u <- par("usr")
  x0 <- u[1] + 0.015 * (u[2] - u[1]); x1 <- u[2] - 0.015 * (u[2] - u[1])
  u_per_in <- abs(diff(grconvertX(c(0, 1), from = "inches", to = "user")))
  slot_u <- min(0.16 * u_per_in, (x1 - x0) / n)
  xc0 <- (x0 + x1) / 2 - n * slot_u / 2
  xs <- xc0 + (seq_len(n) - 0.5) * slot_u
  gap <- slot_u; ymid <- 1.30; v_h <- 0.045
  site_by <- b_sites <- y_sites <- integer(0)
  if (!is.null(ions)) {
    hit <- ions[!is.na(ions$peak), ]
    if (nrow(hit)) {
      hb <- hit[hit$ion == "b", ]
      if (nrow(hb)) b_sites <- sort(unique(as.integer(hb$num)))
      hy <- hit[hit$ion == "y", ]
      if (nrow(hy)) y_sites <- sort(unique(n - as.integer(hy$num)))
      site_by <- sort(unique(c(b_sites, y_sites)))
    }
  }
  op <- par(xpd = TRUE); on.exit(par(op), add = TRUE)
  slot_in <- slot_u / u_per_in
  cexs <- min(8 / 12, max(6 / 12,
           0.9 * slot_in / strwidth("M", cex = 1, units = "inches")))
  cexn <- min(6 / 12, max(5 / 12,
           0.9 * slot_in / strwidth("b88", cex = 1, units = "inches")))
  if (length(site_by)) {
    vx <- (xs[site_by] + xs[site_by + 1L]) / 2
    segments(vx, ymid - v_h, vx, ymid + v_h, col = "black", lwd = 1.33)
  }
  if (length(b_sites)) {
    vxb <- (xs[b_sites] + xs[b_sites + 1L]) / 2
    segments(vxb, ymid - v_h, vxb - gap * 0.5, ymid - v_h, col = "black", lwd = 1.33)
    text(xs[b_sites], ymid - 0.085,
         as.expression(lapply(b_sites, function(k) bquote(b[.(k)]))),
         cex = cexn, col = ION_COL[["b"]], font = 2)
  }
  if (length(y_sites)) {
    vxy <- (xs[y_sites] + xs[y_sites + 1L]) / 2
    segments(vxy, ymid + v_h, vxy + gap * 0.5, ymid + v_h, col = "black", lwd = 1.33)
    text(xs[y_sites + 1L], ymid + 0.085,
         as.expression(lapply(y_sites, function(k) bquote(y[.(n - k)]))),
         cex = cexn, col = ION_COL[["y"]], font = 2)
  }
  if (is.finite(chg_v) && chg_v > 0L)
    text(xc0 - gap * 1.5, ymid, sprintf("%d+ ", chg_v),
         cex = cexs, font = 2, adj = 1, col = "#CA1D52")
  lab_col <- rep("black", n); lab_font <- rep(1, n)
  if (!is.null(mod) && nzchar(mod)) {
    pm <- suppressWarnings(as.integer(sub(",.*$", "",
      strsplit(mod, ";", fixed = TRUE)[[1]])))
    pm <- pm[!is.na(pm) & pm >= 1L & pm <= n]
    lab_col[pm] <- ION_COL[["I"]]; lab_font[pm] <- 2
  }
  text(xs, ymid, strsplit(sq, "")[[1]], cex = cexs, font = lab_font, col = lab_col)
}

draw_single <- function(spec_obj, lines3, chg_v, ann, head_txt) {
  layout(matrix(c(1, 2), 2, 1), heights = c(1.0, 0.35))
  par(oma = c(3, 3, 5.2, 3))
  par(cex = 1)
  mz <- spec_obj$spectrum[[1]]; ints <- spec_obj$spectrum[[2]]
  o <- order(mz, na.last = NA); mz <- mz[o]; ints <- ints[o]
  ints_n <- ints / max(ints[is.finite(ints) & ints > 0])
  xmax <- max(mz, na.rm = TRUE) * 1.03
  # 谱图面板
  par(mar = c(0.5, 4, 3.9, 1))
  plot(0, 0, type = "n", xlim = c(0, xmax), ylim = c(0, 1.15),
       xlab = "", ylab = "", cex.lab = 0.83, cex.axis = 0.83,
       yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i")
  title(ylab = "Relative intensity (%)", cex.lab = 0.83, line = 2.6)
  par(mgp = c(2.0, 0.8, 0))
  axis(2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 100, 20)),
       cex.axis = 0.83, las = 1)
  text(xmax * 0.995, 1.10, ann[[1]], cex = 0.83, adj = c(1, 0.5), xpd = NA)
  text(xmax * 0.995, 1.00, ann[[2]], cex = 0.83, adj = c(1, 0.5), xpd = NA)
  segments(mz, 0, mz, ints_n, lwd = 0.3, col = "grey60")
  th <- all_ions(as.character(cur_seq), as.character(cur_mod), mz, ints)
  if (!is.null(th)) {
    hit <- th[!is.na(th$peak), ]
    if (nrow(hit)) {
      o <- order(ION_PRI[hit$ion], hit$num, hit$chg, na.last = TRUE)
      hit <- hit[o, ]
      spl <- split(hit, hit$peak)
      pk_i <- as.integer(names(spl))
      labs <- vapply(spl, function(d) paste(d$lab, collapse = ""), "")
      cols <- vapply(spl, function(d) ION_COL[[d$ion[1]]], "")
      pmz <- mz[pk_i]; pin <- ints_n[pk_i]
      segments(pmz, 0, pmz, pin, lwd = 0.8, col = cols)
      # 标签不加额外迁移：峰顶上方原位标注
      text(pmz, pin + 0.03, labs, cex = 8/12, srt = 90, col = cols, adj = 0)
    }
  }
  draw_ladder(as.character(cur_seq), th, as.character(cur_mod), chg_v)
  mtext(head_txt, side = 3, line = 4.8, cex = 0.83, adj = 0)
  # ppm 面板
  par(mar = c(2.6, 4, 0.3, 1), xpd = FALSE)
  plot(NA, xlim = c(0, xmax), ylim = c(-20, 20),
       xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i", yaxs = "i")
  title(xlab = expression(italic("m/z")), cex.lab = 0.83, line = 1.8)
  title(ylab = "ppm", cex.lab = 0.83, line = 2.6)
  abline(h = 0, col = "#999999", lwd = 0.6)
  abline(h = c(-20, 20), col = "#BBBBBB", lty = 3, lwd = 0.6)
  if (!is.null(th)) {
    hm <- th[!is.na(th$peak), ]
    if (nrow(hm)) {
      px <- mz[hm$peak]
      py <- (px - hm$mz) / hm$mz * 1e6
      points(px, py, col = ION_COL[hm$ion], pch = 16, cex = 0.5)
    }
  }
  atk <- pretty(c(0, xmax), n = 10)
  atk <- atk[atk >= 0 & atk <= xmax]
  par(mgp = c(2.0, 0.8, 0))
  axis(1, at = atk, cex.axis = 0.83)
  axis(2, at = c(-20, 0, 20), las = 1, cex.axis = 0.83)
  box()
  # 信息行三行全 10 号（信息-序列行零间距；序列-梯、梯-谱正常间距）
  title(lines3[1], outer = TRUE, line = 4.2, cex.main = 0.83)
  title(lines3[2], outer = TRUE, line = 3.1, cex.main = 0.83)
  title(lines3[3], outer = TRUE, line = 2.0, cex.main = 0.83)
}

# ---- 主循环：逐 MGF 拆分 pFind 行，逐谱输出 PDF＋300dpi PNG ---------------------
# pFind-Filtered 按 raw 前缀归组到各 MGF；输出子文件夹＝MGF 文件茎；
# 信息行＝谱图号/打分（含 PSM score）/蛋白质|序列|修饰；16:9 画布
# ---- main loop: split the pFind rows per MGF and write one PDF + 300-dpi PNG
# per spectrum - rows are grouped to their MGF by the raw prefix; outputs go
# into a subfolder named by the MGF stem; info lines = spectrum id / scores
# (including the PSM score) / protein|sequence|modification; 16:9 canvas ----------------------------------------------------------------------
fmt_num <- function(v, dig) {
  x <- suppressWarnings(as.numeric(as.character(v)))
  if (is.na(x)) "-" else sprintf(paste0("%.", dig, "f"), x)
}
n_ok <- 0L
for (fi in seq_along(MGF_FILES)) {
  lu <- read_mgf(MGF_FILES[fi])                           # 逐个 MGF 读取（内存友好）
  if (length(lu$objs) == 0) next
  sub_sel <- sel[which(pf_mgf_idx[sel] == fi)]             # 该 MGF 对应的 pFind 行
  if (!length(sub_sel)) { cat("该 MGF 无选中谱：", MGF_FILES[fi], "\n"); next }
  out_dir <- mgf_stems[fi]                                 # 输出子文件夹
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  cat(sprintf("== MGF %s：%d 谱 ==\n", basename(MGF_FILES[fi]), length(sub_sel)))
  for (i in sub_sel) {

  row <- pf[i, ]
  fn <- as.character(row$File_Name)
  spec_obj <- find_spec(lu, fn)
  if (is.null(spec_obj)) { cat("跳过（MGF 无该谱）:", fn, "\n"); next }
  cur_seq <- as.character(row$Sequence); cur_mod <- as.character(row$Modification)
  mz0 <- spec_obj$spectrum[[1]]; ints0 <- spec_obj$spectrum[[2]]
  chg_v <- tryCatch(as.integer(spec_obj$charge), error = function(e) NA_integer_)
  if (is.na(chg_v)) chg_v <- row$Charge
  if (is.na(chg_v) || chg_v < 1L) chg_v <- 0L
  exp_mh <- tryCatch(as.numeric(spec_obj$pepMass), error = function(e) NA_real_)
  if (is.finite(exp_mh) && chg_v > 1L)
    exp_mh <- exp_mh * chg_v - (chg_v - 1L) * 1.007276
  calc_mh <- NA_real_
  btmp <- try(suppressWarnings(invisible(capture.output(
    bip <- by_ions_position(cur_seq, cur_mod, "Unlabel", aa, aa_tab)))), silent = TRUE)
  if (!inherits(btmp, "try-error") && !is.null(bip) && is.finite(bip$pep_mass))
    calc_mh <- bip$pep_mass
  base_pk <- max(ints0[is.finite(ints0) & ints0 > 0])
  bp_lab <- sub("E\\+0?", "E", sprintf("%.2E", base_pk))
  ann <- c(sprintf("Base peak = %s", bp_lab),
           if (is.finite(exp_mh) && is.finite(calc_mh) && calc_mh > 0)
             sprintf("Mass error = %+.4f Da / %+.2f ppm",
                     exp_mh - calc_mh, (exp_mh - calc_mh) / calc_mh * 1e6)
           else "Mass error = -")
  psm_sc <- "-"
  btmp2 <- try(suppressWarnings(invisible(capture.output(
    b2 <- ions_match(bip, spec_obj$spectrum, 20, "ppm")))), silent = TRUE)
  if (!inherits(btmp2, "try-error") && !is.null(b2) && is.finite(b2$PSM_score_by))
    psm_sc <- sprintf("%.2f", as.numeric(b2$PSM_score_by))
  acc <- "-"; pos_info <- "-"
  if (nzchar(as.character(row$Proteins)))
    acc <- strsplit(as.character(row$Proteins), "/", fixed = TRUE)[[1]][1]
  if (nzchar(as.character(row$Positions))) {
    ft <- strsplit(strsplit(as.character(row$Positions), "/", fixed = TRUE)[[1]][1],
                   ",", fixed = TRUE)[[1]]
    if (length(ft) >= 3) pos_info <- sprintf("%s[%s]%s", ft[2], ft[1], ft[3])
  }
  pre <- sub("\\[.*$", "", pos_info); post <- sub("^.*\\]", "", pos_info)
  head_txt <- sprintf("| %s | (%s)%s(%s) | %s |", acc, pre, cur_seq, post,
                      if (nzchar(cur_mod)) cur_mod else "-")
  lines3 <- c(
    sprintf("Spectrum: %s", fn),
    sprintf("Scores: PSM score = %s | Raw score = %s | Final score = %s",
            psm_sc, fmt_num(row$Raw_Score, 2L), fmt_num(row$Final_Score, 2L)),
    sprintf("Protein: %s | %s | %s", acc, cur_seq,
            if (nzchar(cur_mod)) cur_mod else "-"))
  tag <- if (is.na(row$Scan_No)) i else row$Scan_No
  pdf_f <- file.path(out_dir, sprintf("scan%s.pdf", tag))
  png_f <- file.path(out_dir, sprintf("scan%s.png", tag))
  pdf(pdf_f, width = 12.8, height = 7.2)
  draw_single(spec_obj, lines3, chg_v, ann, head_txt)
  dev.off()
  png(png_f, width = 12.8, height = 7.2, units = "in", res = 300)
  draw_single(spec_obj, lines3, chg_v, ann, head_txt)
  dev.off()
  n_ok <- n_ok + 1L
  if (n_ok %% 20 == 0L) cat(sprintf("进度 %d/%d\n", n_ok, length(sel)))
  }
}
cat(sprintf("完成：输出 %d 张（PDF＋300dpi PNG，子文件夹＝MGF 文件茎）\n", n_ok))
