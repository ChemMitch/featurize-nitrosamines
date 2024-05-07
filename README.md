FeaturizeNitrosamines
=====================

FeaturizeNitrosamines is a java application which computes a set of chemical features and a potency category for a given set of input nitrosamine chemical structures. More information about potency category scores and features can be found here:

Kruhlak, N. L., Schmidt, M., Froetschl, R., Graber, S., Haas, B., Horne, I., Horne, S., King, S. T., Koval, I. A., Kumaran, G., Langenkamp, A., McGovern, T. J., Peryea, T., Sanh, A., Siqueira Ferreira, A., van Aerts, L., Vespa, A., Whomsley, R., 2024. Determining Recommended Acceptable Intake Limits for N-Nitrosamine Impurities in Pharmaceuticals: Development and Application of the Carcinogenic Potency Categorization Approach (CPCA). Regulatory Toxicology and Pharmacology (in press).

https://www.fda.gov/regulatory-information/search-fda-guidance-documents/recommended-acceptable-intake-limits-nitrosamine-drug-substance-related-impurities


Building
========

FeaturizeNitrosamines is a maven project and can be built using the following command:

```
mvn clean package
```

This will produce a packaged executable jar file in the target directory with the filename Featurize-Nitrosamines-<VERSION>.jar.

Running
=======

To run the application, the following command can be used:

```
java -jar target/Featureize-Nitrosamines-0.0.1-SNAPSHOT-jar-with-dependencies.jar -i input_example.txt -o output_example.txt
```

Where example.txt is a tab-delimited list of SMILES strings representing the nitrosamine chemical structures to be evaluated, for example:

```
ID	SMILES	NNO Instance
1	C=CCN(CC=C)N=O	1
2	CN(N=O)C1=CC=NC=C1	1
3	CN(N=O)C1=CC=CN=C1	1
4	OC(=O)C1CCCCN1N=O	1
```

The first column is an ID or name for the structure, which will be repeated in the output file. The second column is a SMILES string for the chemical structure. The third column "NNO Instance" is an optional input column. If present will be carried forward to the output. Otherwise it will be generated.

This will produce an output file input_example.txt of the form:

```
Structure_Name	Nitrosamine Structure	NNO Instance	Potency Category	Potency Score	Alpha-Hydrogens	Alpha-hydrogen score	Tertiary alpha-carbon?	Tertiary alpha-carbon score	Carboxylic acid group anywhere on molecule?	Carboxylic acid group anywhere on molecule score	NNO in pyrrolidine ring?	NNO in pyrrolidine ring score	NNO in 6-membered ring with S?	NNO in 6-membered ring with S score	NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine)?	NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine) score	NNO in morpholine ring?	NNO in morpholine ring score	NNO in a 7-membered ring?	NNO in a 7-membered ring score	Chains of >=5-non-H atoms on both sides of NNO?	Chains of >=5-non-H atoms on both sides of NNO score	EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone)?	EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone) score	EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone)?	EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone) score	Beta-hydroxyl on ONLY one side?	Beta-hydroxyl on ONLY one side score	Beta-hydroxyl on BOTH sides?	Beta-hydroxyl on BOTH sides score	Aryl bonded to alpha-carbon?	Aryl bonded to alpha-carbon score	Methyl group on beta-carbon?	Methyl group on beta-carbon score
1	C=CCNCC=C	1	1	1	2,2	1	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0
2	CNc1ccncc1	1	2	2	0,3	2	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0
3	CNc1cccnc1	1	2	2	0,3	2	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0
4	OC(=O)C1CCCCN1	1	4	8	1,2	3	NO	0	YES	3	NO	0	NO	0	YES	2	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0	NO	0
```

Output Format
=============

The default columns present in the output file are the following:

"Structure_Name"

The supplied first column for an ID/name that was provided in the original input.

"SMILES"

The input structure, parsed and reproduced with potentially added nitrosamine in the case where predictions are used.

"NNO Instance"

Incremental number of nitrosamine group. If provided in the input file, will be repeated, otherwise it will be computed as "1". With the -n flag this will be an incremental number for each new generated nitrosamine.

"Potency Category"

The calculated potency category based on the computed features and the decision tree outlined in the guidance. 

"Potency Score"

The sum of all potency scores based on each feature

Each feature from the program produces a column pair with the first column being the value of the feature (typically YES/NO) and the contributing potency "score" of the calculated feature (typically 0, 1, 2 or 3).

"Alpha-Hydrogens" / "Alpha-hydrogen score"

The count of alpha hydrogens for the nitrosamine group in the form of "<lower>,<higher>", and the potency score for that value.


The other feature columns are:

```
Tertiary alpha-carbon?
Tertiary alpha-carbon score
Carboxylic acid group anywhere on molecule?
Carboxylic acid group anywhere on molecule score
NNO in pyrrolidine ring?
NNO in pyrrolidine ring score
NNO in 6-membered ring with S?
NNO in 6-membered ring with S score
NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine)?
NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered S-containing ring and morpholine) score
NNO in morpholine ring?
NNO in morpholine ring score
NNO in a 7-membered ring?
NNO in a 7-membered ring score
Chains of >=5-non-H atoms on both sides of NNO?
Chains of >=5-non-H atoms on both sides of NNO score
EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone)?
EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, aryl and ketone) score
EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone)?
EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl and ketone) score
Beta-hydroxyl on ONLY one side?
Beta-hydroxyl on ONLY one side score
Beta-hydroxyl on BOTH sides?
Beta-hydroxyl on BOTH sides score
Aryl bonded to alpha-carbon?
Aryl bonded to alpha-carbon score
Methyl group on beta-carbon?
Methyl group on beta-carbon score
```

With definitions consistent with their column names and guidance language. 

Program Options
===============

Running the program with no options will reveal the help menu, and that describes the various input options:

```
usage: java -jar Featureize-Nitrosamines-0.0.1-SNAPSHOT-jar-with-dependencies.jar <options>
 -a,--add-nitrosamines                   add nitrosamines to smiles on
                                         final export at the predicted
                                         site
 -i,--input-file <arg>                   input file (tab-delimited) of
                                         SMILES to process. SMILES should
                                         be second column, 1st column will
                                         be repeated.
 -is,--std-input                         use std input instead of a file
 -nh,--no-headers                        expect no headers on import TSV
 -o,--output-file <arg>                  output file (tab-delimited) of
                                         processed data
 -os,--std-output                        use std output instead of a file
 -p,--print-only                         only print out the input
                                         structures in salt-stripped form
                                         (adding/removing nitrosamines
                                         based on settings) without
                                         generating fingerprint columns
 -pn,--predict-nitrosamines              predict all nitrosamine sites
                                         based on secondary amines and
                                         dimethyl amines, otherwise only
                                         use sites that have atom maps
 -r,--remove-nitrosamines                remove any nitrosamines which are
                                         present on input smiles before
                                         evaluating
 -rh,--remove-headers                    remove headers on export TSV
 -rm,--repress-mapping-site-generation   do NOT generate atom-map
                                         nitrosamine site if none is
                                         present
```



