package gov.fda.gsrs.ndsri;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

public class SaltInfo{
	

	private static String saltList= "MBBZMMPHUWSWHV-BDVNFPICSA-N	14	MEGLUMINE	USE	\r\n" + 
			"GLUUGHFHXGJENI-UHFFFAOYSA-N	1	piperazine	USE	\r\n" + 
			"JUHORIMYRDESRB-UHFFFAOYSA-N	1	BENZATHINE	USE	\r\n" + 
			"FKNQFGJONOIPTF-UHFFFAOYSA-N	302	Sodium	Salt	Cation\r\n" + 
			"BHPQYMZQTOCNFJ-UHFFFAOYSA-N	57	Calcium	Salt	Cation\r\n" + 
			"NPYPAHLBTDXSSS-UHFFFAOYSA-N	57	Potassium	Salt	Cation\r\n" + 
			"JLVVSXFLKOJNIY-UHFFFAOYSA-N	51	Magnesium	Salt	Cation\r\n" + 
			"VTLYFUHAOXGGBS-UHFFFAOYSA-N	35	Iron-3	Salt	Cation\r\n" + 
			"PTFCDOFLOPIGGS-UHFFFAOYSA-N	32	Zinc	Salt	Cation\r\n" + 
			"CWYNVVGOOAEACU-UHFFFAOYSA-N	24	Iron-2	Salt	Cation\r\n" + 
			"RJOJUSXNYCILHH-UHFFFAOYSA-N	11	Gd+3	Salt	Cation\r\n" + 
			"FOIXSVOLVBLSDH-UHFFFAOYSA-N	9	Silver	Salt	Cation\r\n" + 
			"JPVYNHNXODAKFH-UHFFFAOYSA-N	8	Copper	Salt	Cation\r\n" + 
			"LRBQNJMCXXYXIU-PPKXGCFTSA-N	7	TANNATE	Salt	Cation\r\n" + 
			"OEYIOHPDSNJKLS-UHFFFAOYSA-N	7	CHOLINE	Salt	Cation\r\n" + 
			"GBNDTYKAOXLLID-UHFFFAOYSA-N	6	Zirconium	Salt	Cation\r\n" + 
			"JAWGVVJVYSANRY-UHFFFAOYSA-N	6	Cobalt Ion	Salt	Cation\r\n" + 
			"HBBGRARXTFLTSG-UHFFFAOYSA-N	6	Lithium Ion	Salt	Cation\r\n" + 
			"BFGKITSFLPAWGI-UHFFFAOYSA-N	6	Chromium ion	Salt	Cation\r\n" + 
			"XLJKHNWPARRRJB-UHFFFAOYSA-N	5	Cobalt Ion	Salt	cation\r\n" + 
			"BQPIGGFYSBELGY-UHFFFAOYSA-N	5	mercury ion	Salt	cation\r\n" + 
			"JDIBGQFKXXXXPN-UHFFFAOYSA-N	5	bismuth ion	Salt	cation\r\n" + 
			"RJMMFJHMVBOLGY-UHFFFAOYSA-N	2	IN+3	Salt	Cation\r\n" + 
			"CKHJYUSOUQDYEN-YPZZEJLDSA-N	2	Ga+3	Salt	Cation\r\n" + 
			"PSDMOPINLDTFSZ-NJFSPNSNSA-N	1	177Lu+3	Salt	Cation\r\n" + 
			"JPVYNHNXODAKFH-IGMARMGPSA-N	1	Copper 64 ion	Salt	Cation\r\n" + 
			"AFVFQIVMOAPDHO-UHFFFAOYSA-N	32	Mesylate	Salt	A Sulfonic Acid\r\n" + 
			"JOXIMZWYDAKGHI-UHFFFAOYSA-N	6	TOSYLATE	Salt	A Sulfonic Acid\r\n" + 
			"SRSXLGNVWSONIS-UHFFFAOYSA-N	4	BESYLATE	Salt	A Sulfonic Acid\r\n" + 
			"CCIVGXIOQKPBKL-UHFFFAOYSA-N	1	Esylate	Salt	A sulfonic acid\r\n" + 
			"MIOPJNTWMNEORI-GMSGAONNSA-N	1	Camsylate	Salt	A sulfonic acid\r\n" + 
			"KVBGVZZKJNLNJU-UHFFFAOYSA-N	1	NAPSYLIC ACID	Salt	A sulfonic Acid\r\n" + 
			"VEXZGXHMUGYJMC-UHFFFAOYSA-N	314	HCl	Salt	A Strong Acid\r\n" + 
			"QAOWNCQODCNURD-UHFFFAOYSA-N	47	Sulfate	Salt	A Strong Acid\r\n" + 
			"NBIIXXVUZAFLBC-UHFFFAOYSA-N	21	Phosphate	Salt	A Strong Acid\r\n" + 
			"CPELXLSAUQHCOX-UHFFFAOYSA-N	15	HBr	Salt	A Strong Acid\r\n" + 
			"GRYLNZFGIOXLOG-UHFFFAOYSA-N	7	Nitric Acid	Salt	A Strong Acid\r\n" + 
			"QTBSBXVTEAMEQO-UHFFFAOYSA-N	32	Acetate	Salt	A Carboxylic Acid\r\n" + 
			"VZCYOOQTPOCHFL-UPHRSURJSA-N	23	Maleate	Salt	A Carboxylic Acid\r\n" + 
			"FEWJPZIEWOKRBE-JCYAYHJZSA-N	22	Tartrate	Salt	A Carboxylic Acid\r\n" + 
			"VZCYOOQTPOCHFL-OWOJBTEDSA-N	18	Fumeric Acid	Salt	A Carboxylic Acid\r\n" + 
			"KRKNYBCHXYNGOX-UHFFFAOYSA-N	14	CITRATE	Salt	A Carboxylic Acid\r\n" + 
			"KDYFGRWQOYBRFD-UHFFFAOYSA-N	10	SUCCINATE	Salt	A Carboxylic Acid\r\n" + 
			"JVTAAEKCZFNVCJ-UHFFFAOYSA-N	6	Lactic Acid	Salt	A Carboxylic Acid\r\n" + 
			"WLJNZVDCPSBLRP-UHFFFAOYSA-N	4	Pamoate	Salt	A Carboxylic Acid\r\n" + 
			"WPYMKLBDIGXBTP-UHFFFAOYSA-N	4	benzoic acid	salt	A Carboxylic Acid\r\n" + 
			"ODKSFYDXXFIFQN-BYPYZUCNSA-N	3	Arginine	Salt	A Carboxylic Acid\r\n" + 
			"YGSDEFSMJLZEOE-UHFFFAOYSA-N	3	SALICYLIC ACID	Salt	A Carboxylic Acid\r\n" + 
			"RGHNJXZEOKUKBD-SQOUGZDYSA-N	3	Gluconic Acid	Salt	A Carboxylic Acid\r\n" + 
			"CKLJMWTZIZZHCS-REOHCLBHSA-N	3	Aspartic Acid	Salt	A Carboxylic Acid\r\n" + 
			"MUBZPKHOEPUJKR-UHFFFAOYSA-N	3	Oxalic Acid	Salt	A Carboxylic Acid\r\n" + 
			"BJEPYKJPYRNKOW-UHFFFAOYSA-N	2	Malate	Salt	A Carboxylic Acid\r\n" + 
			"BJEPYKJPYRNKOW-REOHCLBHSA-N	2	L-Malate	Salt	A Carboxylic Acid\r\n" + 
			"QIQXTHQIDYTFRH-UHFFFAOYSA-N	2	Stearic Acid	Salt	A Carboxylic Acid\r\n" + 
			"SJJCQDRGABAVBB-UHFFFAOYSA-N	1	XINAFOATE	Salt	A carboxylic acid\r\n" + 
			"DCYGAPKNVCQNOE-UHFFFAOYSA-N	1	TRIFENATATE	Salt	A carboxylic acid\r\n" + 
			"DSLZVSRJTYRBFB-DUHBMQHGSA-N	1	Mucic Acid	Salt	A carboxylic acid\r\n" + 
			"FEWJPZIEWOKRBE-LWMBPPNESA-N	1	S,S-Tartrate	Salt	A carboxylic acid\r\n" + 
			"JYTUSYBCFIZPBE-AMTLMPIISA-N	1	Lactobionate	Salt	A carboxylic acid\r\n" + 
			"IAZDPXIOMUYVGZ-UHFFFAOYSA-N	2	DimethylSulfoxide	Salt	\r\n" + 
			"YBRBMKDOPFTVDT-UHFFFAOYSA-N	1	ERBUMINE	Salt	\r\n" + 
			"XBRDBODLCHKXHI-UHFFFAOYSA-N	1	EPOLAMINE	Salt	\r\n" + 
			"CBTVGIZVANVGBH-UHFFFAOYSA-N	1	AMINOMETHYLPROPANOL	Salt	\r\n" + 
			"XFXPMWWXUTWYJX-UHFFFAOYSA-N	27	Cyano	Exclude All	\r\n" + 
			"RGHNJXZEOKUKBD-SQOUGZDYSA-M	24	Gluconic Anion	Exclude All	\r\n" + 
			"QTBSBXVTEAMEQO-UHFFFAOYSA-M	20	Acetate Anion	Exclude All	\r\n" + 
			"KRKNYBCHXYNGOX-UHFFFAOYSA-K	16	Citrate Anion	Exclude All	\r\n" + 
			"NBIIXXVUZAFLBC-UHFFFAOYSA-K	13	Phosphate Ion	Exclude All	\r\n" + 
			"BVKZGUZCCUSVTD-UHFFFAOYSA-L	13	Acetate Ion	Exclude All	\r\n" + 
			"RWMKKWXZFRMVPB-UHFFFAOYSA-N	12	Silicon ion	Exclude All	\r\n" + 
			"UCKMPCXJQFINFW-UHFFFAOYSA-N	12	Sulfur Anion	Exclude All	\r\n" + 
			"NHNBFGGVMKEFGY-UHFFFAOYSA-N	11	Nitrate Anion	Exclude All	\r\n" + 
			"NBIIXXVUZAFLBC-UHFFFAOYSA-L	10	Phosphate Ion	Exclude All	\r\n" + 
			"LJJFNFYPZOHRHM-UHFFFAOYSA-N	10	Chelator	Exclude All	\r\n" + 
			"DHMQDGOQFOQNFH-UHFFFAOYSA-M	10	Glycinate	Exclude All	\r\n" + 
			"QGZKDVFQNNGYKY-UHFFFAOYSA-O	9	Ammonium	Exclude All	\r\n" + 
			"XPPKVPWEQAFLFU-UHFFFAOYSA-J	8	PYROPHOSPHATE	Exclude All	\r\n" + 
			"ZOKXTWBITQBERF-UHFFFAOYSA-N	8	Molybdenum	Exclude All	\r\n" + 
			"YGSDEFSMJLZEOE-UHFFFAOYSA-M	7	SALICYLATE	Exclude All	\r\n" + 
			"VYPSYNLAJGMNEJ-UHFFFAOYSA-N	6	Silocon Dioxide	Exclude All	\r\n" + 
			"CKLJMWTZIZZHCS-REOHCLBHSA-L	6	Aspartate	Exclude All	\r\n" + 
			"BLRPTPMANUNPDV-UHFFFAOYSA-N	6	Silicon	Exclude All	\r\n" + 
			"XLYOFNOQVPJJNP-UHFFFAOYSA-N	637	Water	Exclude	\r\n" + 
			"AHKZTVQIVOEVFO-UHFFFAOYSA-N	132	Oxygen Ion	Exclude	\r\n" + 
			"VEXZGXHMUGYJMC-UHFFFAOYSA-M	104	Cl-	Exclude	\r\n" + 
			"XLYOFNOQVPJJNP-UHFFFAOYSA-M	93	Hydroxyl	Exclude	\r\n" + 
			"REDXJYDRNCIFBQ-UHFFFAOYSA-N	61	Aluminum	Exclude	\r\n" + 
			"QAOWNCQODCNURD-UHFFFAOYSA-L	39	Sulfate Ion	Exclude	\r\n" + 
			"GPRLSGONYQIRFK-UHFFFAOYSA-N	33	Hydrogen	Exclude	\r\n" + 
			"CPELXLSAUQHCOX-UHFFFAOYSA-M	31	Br-	Exclude	\r\n" + 
			"XMBWDFGMSWQBCA-UHFFFAOYSA-M	14	I-	Exclude	\r\n" + 
			"QGZKDVFQNNGYKY-UHFFFAOYSA-N	7	Ammonia	Exclude	\r\n" + 
			"KRHYYFGTRYWZRS-UHFFFAOYSA-M	6	F-	Exclude	\r\n" + 
			"CZMRCDWAGMRECN-UGDNZRGBSA-N	6	Gluconate	Exclude	\r\n" + 
			"JOXIMZWYDAKGHI-UHFFFAOYSA-M	6	TOSYLATE ANION	Exclude	\r\n" + 
			"SJZRECIVHVDYJC-UHFFFAOYSA-M	6	OXYBATE ION	Exclude	\r\n" + 
			"FZHXIRIBWMQPQF-SLPGGIOYSA-N	5	GLUCOSAMINE	Exclude	\r\n" + 
			"QAOWNCQODCNURD-UHFFFAOYSA-M	3	Sulfate Ion	Exclude	\r\n" + 
			"PIICEJLVQHRZGT-UHFFFAOYSA-N	2	ethylenediamine	Exclude	\r\n" + 
			"LENZDBCJOHFCAS-UHFFFAOYSA-N	7	Tromethamine	Antioxidant	\r\n" + 
			"CIWBSHSKHKDKBQ-JLAZNSOCSA-M	7	Ascorbic Acid	Antioxidant	\r\n" + 
			"KWTSXDURSIMDCE-UHFFFAOYSA-N	7	Amphetamine	AM	\r\n" + 
			"OROGSEYTTFOCAN-DNJOTXNNSA-N	7	Codeine	AM	\r\n" + 
			"MZXXONIZWOBROX-DEOSSOPVSA-N	7	ESOMEPRAZOLE	AM	\r\n" + 
			"BPZSYCZIITTYBL-YJYMSZOUSA-N	6	FORMOTEROL	AM	\r\n" + 
			"ZFXYFBGIUFBOJW-UHFFFAOYSA-N	6	Aminophylline	AM	\r\n" + 
			"RKUNBYITZUJHSG-HEWHGDCOSA-N	6	HYOSCYAMINE	AM	\r\n" + 
			"XUBOMFCQGDBHNK-UHFFFAOYSA-N	5	GatIFLOXACIN	AM	\r\n" + 
			"IUBSYMUCCVWXPE-UHFFFAOYSA-N	4	Metoprolol	AM	";
	
	
	String inchiKey;
	String name;
	String includeType;
	String subType;
	public SaltInfo(String inchiKey, String name, String includeType, String subType){
		this.inchiKey = inchiKey;
		this. name = name;
		this. includeType=includeType;
		this. subType=subType;
	}
	
	public static String getSaltList() {
		return saltList;
	}

	public String getInchiKey() {
		return inchiKey;
	}

	public String getName() {
		return name;
	}

	public String getIncludeType() {
		return includeType;
	}

	public String getSubType() {
		return subType;
	}

	public static Map<String, SaltInfo> defaultSaltMap(){
		Map<String,SaltInfo> saltTable;
	
	
		saltTable=Arrays.stream(saltList.split("\n"))
		      .map(t->t.trim())
		      .map(t->t.split("\t"))
		      .map(a->new SaltInfo(a[0],a[2],a[3],(a.length>=5)?a[4]:""))
		      .collect(Collectors.toConcurrentMap(s->s.inchiKey, s->s));
		return saltTable;
	}
	
}