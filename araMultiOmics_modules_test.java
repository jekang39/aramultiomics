/* 
 * Licensed to the Apache Software Foundation (ASF) under one or more
 *  contributor license agreements.  See the NOTICE file distributed with
 *  this work for additional information regarding copyright ownership.
 *  The ASF licenses this file to You under the Apache License, Version 2.0
 *  (the "License"); you may not use this file except in compliance with
 *  the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 * 
 */
import java.io.*;
import java.util.*;
import java.text.SimpleDateFormat;
import java.util.regex.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.math.BigDecimal;
import org.apache.commons.lang3.Range;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.clustering.*;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.exception.*;
import org.apache.commons.math3.stat.correlation.*;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.util.Combinations;
//javac -cp .:/path/lib/commons-io-2.11.0.jar:/path/lib/commons-lang3-3.12.0.jar:/path/lib/commons-math3-3.6.1.jar araMultiOmics_modules_test.java

//java -cp .:/path/lib/commons-io-2.11.0.jar:/path/lib/commons-lang3-3.12.0.jar:/path/lib/commons-math3-3.6.1.jar araMultiOmics_modules_test>out_araMultiOmics_modules_test.pathway


public class  araMultiOmics_modules_test {

	public static void main(String[] args) throws IOException,InterruptedException,Exception{

		//instantiate main
		araMultiOmics_modules_test aa = new araMultiOmics_modules_test();

		//instantiate inner class dni
	    	//araMultiOmics_modules_test.dni_module dni_chip_atac_seq = new araMultiOmics_modules_test.dni_module();

		//instantiate inner class dna_tf
	    	//araMultiOmics_modules_test.dna_tf_module dna_tf = new araMultiOmics_modules_test.dna_tf_module();
	    	
	    	//instantiate inner class DNA_struc_metal_binding  
		//araMultiOmics_modules_test.DNA_struc_metal_binding dna_struc_metal_bind = new araMultiOmics_modules_test.DNA_struc_metal_binding();

	      //instantiate inner class ddi interface  
		//araMultiOmics_modules_test.ddi_interface_anno_module ddi_anno = new araMultiOmics_modules_test.ddi_interface_anno_module();
	    
		//instantiate inner class pathway analysis 
		//araMultiOmics_modules_test.pathway_anal_module path_analysis = new araMultiOmics_modules_test.pathway_anal_module();

		//instantiate inner class Muti-omics association
		araMultiOmics_modules_test.multiomics_asso_module multi_asso = new araMultiOmics_modules_test.multiomics_asso_module();
		
		//instantiate inner class PO_GO
		//araMultiOmics_modules_test.PO_GO po_slim_module = new araMultiOmics_modules_test.PO_GO();
		
		//instantiate inner medicago_cog_module
		//araMultiOmics_modules_test.medicago_cog_module medi_cog = new araMultiOmics_modules_test.medicago_cog_module();


		////inner class functions

		//dni_chip_atac_seq.get_how_many_peak();
		//dni_chip_atac_seq.get_how_many_peak_spec();
		//dni_chip_atac_seq.init_ocr_anal();
		//dni_chip_atac_seq.map_gene_pcsd_region_state();
		//dni_chip_atac_seq.map_sra_exp();
		//dni_chip_atac_seq.map_gene_sra_exp();
		//dni_chip_atac_seq.map_sra_to_tissue();
		//dni_chip_atac_seq.map_gene_to_tissue();
		
		//dni_chip_atac_seq.map_motif_acc_to_domain();
		//dni_chip_atac_seq.calc_atac_width_phloem();
		//dni_chip_atac_seq.create_phloem_atac_seq();
		
		//dna_tf.map_atxg_to_gene_name();
	        //dna_tf.map_atxg_to_tf_name();
	        //dna_tf.link_tf_to_fam();
	        
	        //dna_struc_metal_bind.get_ara_gquad_ortho();
		//dna_struc_metal_bind.map_gene_to_k_dependent_qquad();
		//dna_struc_metal_bind.map_gene_to_k_PDS_dependent_qquad();
		//dna_struc_metal_bind.map_medi_string_uniprot_to_gff_transcript();
		//dna_struc_metal_bind.map_ara_string_uniprot_to_gff_transcript();
		//dna_struc_metal_bind.map_ara_string_uniprot_to_gff_id();
		//dna_struc_metal_bind.map_ara_atxg_to_dna_g4();
		//dna_struc_metal_bind.map_medi_mtr_to_dna_g4();
		//dna_struc_metal_bind.map_ara_atxg_to_gquad_prog();
		//dna_struc_metal_bind.map_medi_mtr_to_gquad_prog();
		//dna_struc_metal_bind.map_ara_atxg_to_prot_metal_binding();
		//dna_struc_metal_bind.map_medi_mtr_to_prot_metal_binding();
	        
		//ddi_anno.map_domain_pfam_3did();
    		//ddi_anno.map_cat_stat();
    		//ddi_anno.map_pdbid_to_pfam();
    		//ddi_anno.map_pdb_to_uniprot();
    		//ddi_anno.map_pfam_to_tair_gene();
    		//ddi_anno.map_atxg_to_uniprot();
    		//ddi_anno.make_string_db_tair();
    		//ddi_anno.make_tf_cofactor_list_tair();
   		
    		
		path_analysis.retrieve_plantcyc_comp_type();
      		path_analysis.map_gene_name_to_pmn_acc();
		path_analysis.map_product_to_enzrxn_id();
		path_analysis.map_compound_id_to_features();
		path_analysis.map_enzrxn_id_to_features();
      		path_analysis.from_pwy_to_acc();
      		path_analysis.map_rxn_id_to_enzrxn_id();
      		
		//multi_asso.map_pheno_to_path_study16();
		//multi_asso.map_pheno_to_path();
		//multi_asso.get_isoform_flag();
		//multi_asso.map_gene_to_pheno_twas();
		
		po_slim_module.get_cluster_genes_po();
		po_slim_module.get_cluster_genes_slim();
		//po_slim_module.map_go_term_type_to_term();
		//po_slim_module.map_gene_to_ipd3();
		
 		//medi_cog.map_medi_string_uniprot_to_gff_transcript();
  		//medi_cog.map_ara_string_uniprot_to_gff_transcript();
  		//medi_cog.map_ara_string_uniprot_to_gff_id();
 		//medi_cog.map_ara_to_cog();
  		//medi_cog.map_medi_to_cog();
  		//medi_cog.get_cog_ara_medi_intersection();
  		//medi_cog.match_cogs_between_ara_medi_intersection();
  		//medi_cog.get_ara_uniprot_module_e_cog();
  		//medi_cog.get_module_e_medi_cog();
  		//medi_cog.get_medi_deg();
  		//medi_cog.map_medi_to_inter_part();

		//main class functions
		
		//test dni module
		
		//aa.map_gene_to_tissue_specific_enhancer();
		//aa.test_dni();
		//aa.test_init_ocr();
		
		//test ddi module
		//aa.test_ddi();
		
		//test pathway
		aa.test_compound();
		
}//main function

//tests start

String[] homologue_path = {"HAI","CYP76C","PAL","CSLG","VTE","GUS"};
String[][] homologue_genes = {{"HAI1","HAI2","HAI3"},{"CYP76C1","CYP76C2","CYP76C3","CYP76C4","CYP76C5","CYP76C6","CYP76C7","CYP76C8P"},{"PAL1","PAL2","PAL3","PAL4"},{"CSLG1","CSLG2","CSLG3"},{"VTE1","VTE5","VTE6","VTE7"},{"GUS1","GUS2","GUS3"}};

String[][] homologue_atxg = {{"AT5G59220","AT1G07430","AT2G29380"},{"AT2G45560","AT2G45570","AT2G45580","AT2G45550","AT1G33730","AT1G33720","AT3G61040","AT3G60955"},{"AT2G37040","AT3G53260","AT5G04230","AT3G10340"},{"AT4G24010","AT4G24000","AT4G23990"},{"AT4G32770","AT5G04490","AT1G78620","AT5G39220"},{"AT5G61250","AT5G07830","AT5G34940"}};

String[][] homologue_module = {{"F","B","B"},{"R","E","NA","R","NA","NA","NA","NA"},{"Q","D","NA","B"},{"D","E","NA"},{"E","B","R","NA"},{"D-A-B","E","Q-B"}};

///
String[] map_go = {
"ecisbinding:enables transcription cis-regulatory region binding","imetabolic:involved in cinnamic acid biosynthetic process","eprotbinding:enables protein binding","eenzyme:enables phenylalanine ammonia-lyase activity","awaterstress:acts upstream of or within response to water deprivation","asaltstress:acts upstream of or within response to salt stress","aresplipid:acts upstream of or within response to lipid","aresplight:acts upstream of or within response to light stimulus","aunigrowth:acts upstream of or within regulation of unidimensional cell growth","areproduct:acts upstream of or within regulation of reproductive process","amultiorgaprocess:acts upstream of or within regulation of multicellular organismal process","evoc:acts upstream of or within phenylpropanoid metabolic process","evoc:acts upstream of or within phenylpropanoid biosynthetic process","aorgcatabol:acts upstream of or within organic substance catabolic process","alipidmeta:acts upstream of or within lipid biosynthetic process","adefen:acts upstream of or within defense response to bacterium","adefen:acts upstream of or within defense response","aorgcatabol:acts upstream of or within cellular catabolic process","adifferentiation:acts upstream of or within cell differentiation","loccyto:located in cytoplasm ","evoc:involved in L-phenylalanine catabolic process","astress:acts upstream of or within response to wounding","astress:acts upstream of or within response to oxidative stress","evoc:acts upstream of or within L-phenylalanine catabolic process","areproduct:acts upstream of or within gametophyte development","aorgcatabol:acts upstream of or within salicylic acid catabolic process","imetabolic:acts upstream of or within small molecule metabolic process","aresplight:acts upstream of or within response to UV-B","aorgcatabol:acts upstream of or within lignin catabolic process","awaterstress:acts upstream of or within drought recovery","locnucl:is active in nucleus ","idephos:involved in peptidyl-threonine dephosphorylation","idephos:enables protein serine/threonine phosphatase activity","idephos:enables myosin phosphatase activity","adifferentiation:acts upstream of or within release of seed from dormancy","adifferentiation:acts upstream of or within positive regulation of seed germination","aposga:acts upstream of or within positive regulation of gibberellic acid mediated signaling pathway","adifferentiation:acts upstream of or within negative regulation of seed dormancy process","anegaba:acts upstream of or within negative regulation of abscisic acid-activated signaling pathway","locchlo:located in chloroplast","locnucl:located in nucleus","loccyto:located in cytosol","locgolgimem:located in cis-Golgi network membrane","idephos:involved in protein dephosphorylation","astoma:acts upstream of or within stomatal movement","aposaba:acts upstream of or within response to abscisic acid","aleafsenesc:acts upstream of or within leaf senescence","achloro:acts upstream of or within chloroplast organization","locgolgi:located in Golgi apparatus","aorga:acts upstream of or within organelle organization","roxygen:enables oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","roxygen:enables monooxygenase activity","emetalbind:enables iron ion binding ","emetalbind:enables heme binding","areproduct:acts upstream of or within pollen development"
};
////
//
String[] map_atxg_go_var = {"AT1G01360  eprotbinding","AT1G07430  adifferentiation","AT1G07430  anegaba","AT1G07430  aposga","AT1G07430  eprotbinding","AT1G07430  idephos","AT1G07430 locchlo","AT1G07430  locnucl","AT1G55270  eprotbinding","AT1G73000  eprotbinding","AT2G26040  eprotbinding","AT2G29380  eprotbinding","AT2G29380  idephos","AT2G29380 locnucl","AT2G29380  locnucl","AT2G35900  aorga","AT2G37040  adefen","AT2G37040  aorgcatabol","AT2G37040  areproduct","AT2G37040  aresplight","AT2G37040  astress","AT2G37040  awaterstress","AT2G37040  eenzyme","AT2G37040  eprotbinding","AT2G37040  evoc","AT2G37040  imetabolic","AT2G37040 loccyto","AT2G38310  eprotbinding","AT2G40330  eprotbinding","AT2G45560  aresplight","AT2G45560  emetalbind","AT2G45560  emetalbind","AT2G45560 locchlo","AT2G45560  roxygen","AT3G10340  adefen","AT3G10340  adefen","AT3G10340  adifferentiation","AT3G10340  alipidmeta","AT3G10340  amultiorgaprocess","AT3G10340  aorgcatabol","AT3G10340  aorgcatabol","AT3G10340  areproduct","AT3G10340  aresplight","AT3G10340  aresplipid","AT3G10340  asaltstress","AT3G10340  aunigrowth","AT3G10340  awaterstress","AT3G10340  eprotbinding","AT3G10340  evoc","AT3G10340  imetabolic","AT3G10340 loccyto","AT3G15180  areproduct","AT3G15180  awaterstress","AT3G53260  adefen","AT3G53260  astress","AT3G53260  eenzyme","AT3G53260  eprotbinding","AT3G53260  evoc","AT3G53260  imetabolic","AT3G53260 loccyto","AT3G54810  ecisbinding","AT4G01026  eprotbinding","AT4G17490  ecisbinding","AT4G17870  eprotbinding","AT4G18620  eprotbinding","AT4G23980  ecisbinding","AT4G27920  eprotbinding","AT4G36620  ecisbinding","AT5G02840  ecisbinding","AT5G05440  eprotbinding","AT5G11260  ecisbinding","AT5G15210  ecisbinding","AT5G15600  imetabolic","AT5G44210  ecisbinding","AT5G45860  eprotbinding","AT5G45870  eprotbinding","AT5G46790  eprotbinding","AT5G47220  ecisbinding","AT5G53160  eprotbinding","AT5G59220  achloro","AT5G59220  aleafsenesc","AT5G59220  anegaba","AT5G59220  aposaba","AT5G59220  astoma","AT5G59220  awaterstress","AT5G59220  eprotbinding","AT5G59220  idephos","AT5G59220 locchlo","AT5G59220  loccyto","AT5G59220 locgolgi","AT5G59220  locgolgimem","AT5G59220  locnucl","AT5G59220  locnucl","AT5G60690  ecisbinding"};

LinkedHashMap<String, String> go_var_to_desc;
void map_go_var(){

	go_var_to_desc = new LinkedHashMap();

	Pattern colon_pattern = Pattern.compile(":");
	
	for(int z = 0; z < map_go.length; z++){
	
		String temp = map_go[z].trim();
		String var = temp.substring(0, temp.indexOf(":"));
		String desc = temp.substring(temp.indexOf(":") +1,temp.length());
		go_var_to_desc.put(var, desc);
	
	}
	
}

LinkedHashMap<String, LinkedList<String>> atxg_to_go_var;
void map_atxg_to_go_var(){

	atxg_to_go_var = new LinkedHashMap();

	Pattern space_pattern = Pattern.compile("\\s+");
	
	for(int z = 0; z < map_atxg_go_var.length; z++){
	
		String temp = map_atxg_go_var[z].trim();
		String[] split_var = space_pattern.split(temp);
		String atxg = split_var[0].trim();
		String var = split_var[1].trim();
		if(atxg_to_go_var.containsKey(atxg)){

			LinkedList<String> exist =atxg_to_go_var.get(atxg);

			if(!exist.contains(var)){
				exist.add(var);
			}
		
								
			atxg_to_go_var.put(atxg,exist);
		}else{

			LinkedList<String> exist = new LinkedList();
			exist.add(var);
			atxg_to_go_var.put(atxg,exist);
								
		}
	
	}
	
}


void header() throws IOException{
 	String header = "";
	for(int j = 0; j < dni_module.gene_region.length; j++){

		for(int e = 0; e < dni_module.chip_tissue.length; e++){
						
			header = header + dni_module.gene_region[j] + "_" + dni_module.chip_tissue[e] + ",";
		}
	}

	System.out.println("enhancer specificity header:"+header);
}

LinkedHashMap<String,Boolean[][][]> gene_to_tissue_specific_enhancer;

void map_gene_to_tissue_specific_enhancer() throws IOException {
	System.out.println("map_gene_to_tissue_specific_enhancer():" );

	gene_to_tissue_specific_enhancer = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern semi_colon_pattern = Pattern.compile(";");
	Pattern digit_pattern = Pattern.compile("\\d");

	String path = "data/ewas.tab.2.study16.2.sort.homolog.ann.short.inter.chip_hub.csv";

	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
						
		/*
3	19744454	19744455	.	.	.	CG	Mo98	molybdenum concentration	AT3G53260	3	19744320	19747804	AT3G53260	Exon	Enhancer	SRX7780581;SRX7780580;SRX7780573;SRX7780572;SRX7780559;SRX7780574;SRX7780560;SRX7780566;SRX7780555;SRX7780567;SRX7780579;SRX7780554;SRX7780553;SRX7780565;SRX6784998;SRX6784997;SRX7780567;SRX7780580;SRX7780573;SRX7780554;SRX7780579;SRX7780574;SRX7780555;SRX7780565;SRX7780553;SRX7780559;SRX7780560;SRX4101263;SRX7780574;SRX8861315;SRX3348143;SRX6785004;SRX9132138;SRX8840342;SRX8861314;SRX3006644;SRX3348146;SRX3348141;SRX8861307;SRX6785002;SRX6785008;SRX3006646;SRX1204335;SRX8843252;SRX8840343;SRX8861325;SRX8861305;SRX9132136;SRX3348150;SRX3348147;SRX8861331;SRX8843253;SRX9132137;SRX3348152;SRX3348138;SRX6784996;SRX3348151;SRX3348144;SRX391972;SRX3348149;SRX10771698;SRX3348145;SRX3348148;SRX3348137;SRX391964;SRX8861330;SRX6785007;SRX3041699;SRX2000803;SRX6081144;SRX4382141;SRX1204325;SRX111008;SRX3348142;SRX6785020;SRX8861324;SRX6785001;SRX7442797;SRX6785003;SRX8861329;SRX391962;SRX8861322;SRX8861332;SRX9770829;SRX8861297;SRX8861323;SRX6785006;SRX8861309;SRX111005;SRX111011;SRX9770831;SRX6785005;SRX7442794;SRX3348140;SRX391959;SRX3006645;SRX3041698;SRX4382140;SRX8861333;SRX111007;SRX10089717;SRX9770823;SRX111004;SRX9770864;SRX9770821;SRX111006;SRX8861313;SRX10089718;SRX9770822;SRX4382142;SRX8843251;SRX9770856;SRX3348139;SRX10089719;SRX391963;SRX9770824;SRX9770784;SRX8861304;SRX10089722;SRX3041697;SRX4382143;SRX9770855;SRX9074821;SRX2000801;SRX7442795;SRX9770853;SRX391966;SRX2000808;SRX9770863;SRX5052459;SRX3503848;SRX4382145;SRX5052462;SRX3503850;SRX9770799;SRX10089720;SRX111010;SRX7442793;SRX8861319;SRX5036315;SRX1096550;SRX7442791;SRX8861318;SRX391960;SRX2000804;SRX3503847;SRX5052461;SRX2000799;SRX5534551;SRX5052457;SRX2528910;SRX6081143;SRX8861288;SRX4382144;SRX8861291;SRX8861295;SRX5052463;SRX7442790;SRX5088431;SRX5036313;SRX9770868;SRX8861326;SRX8861301;SRX9770792;SRX2000805;SRX8861303;SRX9770818;SRX9770830;SRX9770854;SRX8861327;SRX8861299;SRX8861334;SRX8861328;SRX9770866;SRX8861306;SRX9770865;SRX9770791;SRX9770789;SRX8861298;SRX9770867;SRX2000806;SRX9770832;SRX8861316;SRX2000807;SRX8861317;SRX9770801;SRX9770775;SRX8861329;SRX7780573;SRX7780572;SRX7780580;SRX7780581;SRX7780560;SRX895361;SRX7780574;SRX7780579;ERX3444803;ERX3444764;SRX391972;SRX895271;SRX8840342;SRX8861325;ERX3444765;ERX3444760;ERX3444767;ERX3444799;SRX8840343;SRX2528906;SRX391994;SRX10089717;SRX1204325;SRX3040869;SRX3348151;SRX2000808;SRX3348141;SRX10771698;SRX3348144;SRX895327;SRX895322;SRX2528910;SRX2000803;ERX3444761;SRX895321;SRX3348145;ERX3444804;SRX895364;SRX895336;SRX391964;SRX895332;SRX7780566;SRX895266;SRX8861305;SRX8861330;SRX2528909;SRX9770862;SRX9074821;SRX895323;SRX9770854;SRX5413160;SRX3348148;SRX5413169;SRX277579;SRX3348147;SRX3040867;SRX3006644;SRX391995;SRX10089719;SRX8861329;SRX3041701;SRX5036315;SRX391968;SRX5413168;SRX3348139;SRX3348146;SRX895366;SRX3348150;SRX4382145;SRX3040866;SRX391971;SRX9132138;SRX111011;SRX111008;ERX3444762;SRX9770863;SRX2528905;SRX895326;SRX391996;SRX111010;SRX4382143;SRX5413161;SRX111004;SRX3006646;SRX9770824;SRX111005;SRX2000800;SRX3348152;SRX3348143;SRX3348138;SRX8861302;SRX1204335;SRX8861324;SRX8861327;SRX3503849;SRX8861313;SRX895330;SRX895237;SRX895220;SRX6785007;SRX391991;SRX8861315;SRX4382142;SRX3041702;SRX10089722;SRX5534550;SRX6785006;SRX8861326;SRX5036313;SRX2528907;SRX3041700;SRX10089718;SRX10089720;SRX1096549;SRX277582;ERX3444806;SRX2528908;SRX3041698;SRX9770821;SRX3006645;SRX3006647;SRX9770796;SRX8861300;SRX895270;SRX895316;SRX391993;SRX7442790;SRX2000807;SRX3348142;SRX2528907;SRX10089723;SRX3348137;SRX895328;SRX6081144;SRX7442794;SRX3348149;SRX7442791;SRX9770793;SRX277581;SRX4382144;SRX895331;SRX9770791;SRX8861314;SRX8861290;SRX9132137;SRX8861328;SRX111007;SRX391967;SRX5036314;SRX4916679;SRX7442792;SRX9770864;SRX6785002;SRX2311144;SRX8861310;SRX1098135;SRX6785004;SRX6785003;SRX895324;SRX8861332;SRX3503850;SRX9132136;SRX277580;SRX3348140;SRX7442797;SRX2000804;SRX894596;SRX2000801;SRX895317;SRX4916680;SRX9770794;SRX9770814;SRX7442795;SRX3041699;SRX4382141;SRX4382140;SRX4916684;SRX391970;SRX8861319;SRX1098136;SRX1096548;SRX111009;SRX9074820;SRX8861311;SRX895244;SRX3041697;SRX895329;SRX391988;SRX9770822;SRX10089721;SRX9074822;SRX8861292;ERX3444802;SRX895359;SRX2000799;SRX5534551;SRX895334;SRX7442793;SRX3503847;SRX9770825;SRX2000802;SRX6785008;SRX6785001;SRX9770792;SRX2000805;SRX8861303;SRX391962;SRX8861299;SRX895269;SRX9770787;SRX9770789;SRX7780565;SRX9770797;SRX9770812;SRX9770806;SRX8861298;ERX3444807;SRX895318;SRX8861331;SRX9770837;SRX8861316;SRX8861322;SRX895319;SRX2000806;SRX9770823;SRX5088431;SRX3040870;SRX9770811;SRX9770790;SRX8861296;SRX2528907;SRX895212;ERX3444801;SRX9770840;SRX5088433;SRX9770829;SRX9770815;SRX8861295;SRX8861301;ERX3444808;SRX8861333;SRX9770838;SRX9770800;SRX8861323;SRX9770817;SRX9770832;SRX9770831;SRX8861288;SRX8861289;ERX3444770;SRX8861297;SRX8861308;SRX391989;SRX9770839;SRX10089724;SRX4916683;SRX4916678;ERX3444805;SRX391990;SRX4916677;SRX8861312;SRX3503848;SRX8861291;SRX8861307;SRX391963;SRX6784996;ERX3444769;SRX9770819;SRX8861334;SRX5052461;SRX8861318;SRX391987;SRX1096550;SRX9770813;SRX7442796;SRX9770782;SRX9770833;ERX3444768;SRX391997;SRX9770865;SRX6081143;SRX391986;SRX9770795;SRX111006;SRX3040868;SRX9770805;SRX9770809;SRX9770786;SRX8861304;SRX5088432;SRX9770826;SRX391966;SRX9770827;SRX6785005;SRX8861293;SRX9770855;SRX391960;SRX391992;SRX2000812;SRX8861317;SRX9770807;SRX9770818;SRX391959;SRX1098134;SRX2528906;SRX895355;SRX895226;SRX9770820;SRX391969;SRX1098137;SRX5088430;SRX2311146;SRX5052457;SRX7780581;SRX9770845;ERX3444797;ERX3444800;SRX391965;SRX7780555;SRX9770834;ERX3444798;SRX7780554;SRX9770867;SRX9770816;SRX7780580;SRX1096551;SRX9770856;SRX6785020;SRX9770794;SRX895333;SRX111009;SRX111004;SRX9770853;SRX9770829;SRX1098138;SRX9770791;SRX9770837;SRX9770866;SRX8861321;SRX9770789;SRX9770868;SRX9770839;SRX7780574;SRX9770796;SRX9770830;SRX9770809;SRX391961;SRX8861301;SRX9770828;SRX9770831;SRX9770795;SRX8861292;SRX9770840;SRX9770838;SRX9770793;SRX8861317;SRX9770819;SRX9770785;SRX9770832;SRX9770817;SRX9770805;SRX9770812;SRX5088433;SRX9770827;SRX9770810;SRX8861288;SRX9770780;SRX9770815;SRX9770811;SRX9770800;SRX9770777;SRX9770788;SRX9770786;SRX9770846;SRX9770804;SRX9770784;SRX8844363;SRX391997;SRX9770790;SRX8861295;SRX9770855;SRX9770847;SRX9770845;SRX9770773;SRX9770783;SRX9770799;SRX895334;SRX9770835;SRX9770836;SRX5052462;SRX9770818;SRX8844365;SRX9770781;SRX2000812;SRX9770797;ERX3444804;ERX3444802;ERX3444798;ERX3444805;SRX895324;ERX3444800;ERX3444807;ERX3444808	3	19745694	19750228	2110	1

==

2	18779113	18779114	.	.	.	CG	Co59	cobalt concentration	ID=AT2G45560;Note=protein_coding_gene;Name=AT2G45560
2	18779113	18779114	.	.	.	CG	K39	potassium concentration	ID=AT2G45560;Note=protein_coding_gene;Name=AT2G45560


*/

		String line = lines.get(q).trim();
		String[] tab_split = tab_pattern.split(line);
		String cur_region = "";
		
		for(int t = 0; t < dni_module.gene_region.length; t++){
		
			if(line.contains(dni_module.gene_region[t])){
			
				cur_region = dni_module.gene_region[t];
			}
		}
		
		String cur_type = "";
		
		
		for(int t = 0; t < dni_module.type.length; t++){
		
			if(line.contains(dni_module.type[t])){
			
				cur_type = dni_module.type[t];
				
			}
		}
		String sra = tab_split[16].trim();
		String cur_gene = tab_split[9].trim();
		String cur_start = tab_split[1].trim();
		//System.out.println("error start_list.size():"+cur_gene+":"+line);
		LinkedList<String> start_list = dni_module.peak_start.get(cur_gene);
		//System.out.println("start_list.size():"+cur_gene + ":" + start_list.size());
		
		if(start_list != null){
		
			Boolean[][][] val = new Boolean[dni_module.gene_region.length][][];

			for(int j = 0; j < val.length; j++){

				val[j] = new Boolean[start_list.size()][];
				
				for(int l = 0; l < val[j].length; l++){
					val[j][l] = new Boolean[dni_module.chip_sra_to_tissue.length];

					for(int e = 0; e < dni_module.chip_sra_to_tissue.length; e++){
						
						val[j][l][e] = new Boolean(false);

						
					}
				}
			}
			//System.out.println("error 1674:"+val.length+":"+ dni_module.chip_sra_to_tissue.length);
			String header = "";
			for(int j = 0; j < dni_module.gene_region.length; j++){

				for(int e = 0; e < dni_module.chip_sra_to_tissue.length; e++){
						
					header = header + dni_module.gene_region[j] + "_" + dni_module.chip_tissue[e] + ",";
				}
			}

			System.out.println("enhancer specificity header:"+header);
			
			int region_int = -1;

			for(int i = 0; i < dni_module.gene_region.length; i++){
			
				if(cur_region.equals(dni_module.gene_region[i])){

					region_int = i;
					break;
				}
			}
			
			int start_int = -1;
			for(int e = 0; e < start_list.size(); e++){
				if(cur_start.equals(start_list.get(e))){

					start_int = e;
					break;

				}
			}
			//System.out.println("error 1706:"+region_int+":"+ start_int);
			//System.out.println("cur_gene:" + cur_gene + ":" + cur_start );
			
			if(region_int != -1 && start_int != -1){

				String[] sra_split = semi_colon_pattern.split(sra);

				for(int e = 0; e < sra_split.length; e++){

					String cur_sra = sra_split[e].trim();
					System.out.print("cur_sra:"+cur_sra);
					for(int m = 0; m < dni_module.chip_sra_to_tissue.length; m++){

						//System.out.println("error 1706:"+region_int+":"+ start_int);
						if(dni_module.chip_sra_to_tissue[m] != null){

							for(int n = 0; n < dni_module.chip_sra_to_tissue[m].length; n++){
				
								if(cur_sra.equals(dni_module.chip_sra_to_tissue[m][n])){

									val[region_int][start_int][m] = new Boolean(true);

									break;
								}
							}
							
							//System.out.println("error 1729:"+region_int+":"+ start_int+":"+m);
							if(val[region_int][start_int][m] != null){
							
							
								if(val[region_int][start_int][m].booleanValue()){

									break;
								}
							}
						}
						
					}//m
				}//e


				
				if(gene_to_tissue_specific_enhancer.containsKey(cur_gene)){
						
					Boolean[][][] exist = gene_to_tissue_specific_enhancer.get(cur_gene);

					for(int x = 0; x < val[region_int][start_int].length; x++){
					
						if(val[region_int][start_int][x].booleanValue()){
							exist[region_int][start_int][x] = new Boolean(true);
										
								
							
						}
					}
						
					gene_to_tissue_specific_enhancer.put(cur_gene,exist);
							
				}else{
					gene_to_tissue_specific_enhancer.put(cur_gene,val);
				}
			}//not -1
		}//not null
		else{
		
			System.out.println("peak null error:" + line);
			System.out.println("peak null error:" + cur_gene + ":" + cur_start);
		}

	}//lines
	
	List<String> keys = new LinkedList<String>(gene_to_tissue_specific_enhancer.keySet());
	for(int i = 0; i < keys.size(); i++){

		String acc = keys.get(i);
		LinkedList<String> start_list = dni_module.peak_start.get(acc);
		Boolean[][][] temp_tissue = gene_to_tissue_specific_enhancer.get(acc);
		for(int q = 0; q < dni_module.gene_region.length; q++){
			
			for(int w = 0; w < start_list.size(); w++){
					
				for(int r = 0; r < dni_module.chip_tissue.length; r++){
					
					System.out.println("gene_to_tissue_specific_enhancer:" +  acc + ":" + dni_module.gene_region[q] + ":" + start_list.get(w) + ":" + dni_module.chip_tissue[r] + ":" + temp_tissue[q][w][r].booleanValue());

				}
			}
		}
		
	}

}//method	


void test_dni() throws IOException{

	System.out.println("test:" );

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");

	//the number of the compound type was 1881

	

	for(int w = 0; w < homologue_atxg.length; w++){
		
		for(int q= 0; q < homologue_atxg[w].length; q++){


				String gene = homologue_atxg[w][q];
				Integer[][] pcsd = dni_module.gene_to_region_state.get(gene);
				
				String pcsd_str = "";
				
				if(pcsd != null){
				
					for(int r = 0; r < pcsd.length; r++){
					
						for(int t = 0; t < pcsd[r].length; t++){
						
							pcsd_str = pcsd_str + pcsd[r][t].intValue() + ",";
						}
					}
				}else{
				
					for(int r = 0; r < dni_module.gene_region.length; r++){
					
						for(int t = 0; t < dni_module.state.length; t++){
						
							pcsd_str = pcsd_str +  "NA,";
						}
					}
				
				}
				
				String width = "";
				Boolean[] val = dni_module.atac_to_width.get(gene);
				
				if(val == null){
				
					
					for(int r = 0; r < dni_module.ranges.length; r++){
					
						width = width + "null,";
					}
					
				}else{
				
					for(int r = 0; r < val.length; r++){
					
						width = width + val[r].booleanValue() + ",";
					}
				}
				
				System.out.println("width input:" + gene + "," + homologue_module[w][q] + "," + pcsd_str + width + homologue_path[w] );
				//String ewas = multi_asso_test(gene);
				
				
				
				Boolean[][][][] gene_tissue = dni_module.gene_to_tissue.get(gene);
				String tissue_str = "";
				if(gene_tissue == null){
				
					for(int j = 0; j < dni_module.type.length; j++){//type

						for(int l = 0; l < dni_module.gene_region.length; l++){//region
							
							for(int e = 0; e < dni_module.chip_tissue.length; e++){

								tissue_str = tissue_str + "NA,";
								

							}
						}
					}
				}else{
				
					tissue_str = "";
					boolean[][][] cur_gene_tissue = new boolean[dni_module.type.length][dni_module.gene_region.length][dni_module.chip_tissue.length];
					
					for(int j = 0; j < gene_tissue.length; j++){//type

						for(int l = 0; l < gene_tissue[j].length; l++){//region
							
							for(int k = 0; k < gene_tissue[j][l].length; k++){//start
							
								for(int e = 0; e < dni_module.chip_tissue.length; e++){

									if(gene_tissue[j][l][k][e].booleanValue()){
										cur_gene_tissue[j][k][e] = true;
									}
								}

							}
						}
					}
					
					
					
					for(int j = 0; j < gene_tissue.length; j++){//type

						for(int l = 0; l < gene_tissue[j].length; l++){//region
							
							for(int e = 0; e < dni_module.chip_tissue.length; e++){

								tissue_str = tissue_str + cur_gene_tissue[j][l][e]  + ",";	

							}
						}
					}
				}
				
				System.out.println("tissue input:" + gene + "," + homologue_module[w][q]+ "," + pcsd_str + tissue_str + homologue_path[w] );
				
				Boolean[][][] enhancer_spec = gene_to_tissue_specific_enhancer.get(gene);
				
				String enhancer_str = "";
				
				if(enhancer_spec == null){
				
					for(int j = 0; j < dni_module.gene_region.length; j++){

						for(int e = 0; e < dni_module.chip_sra_to_tissue.length; e++){
								
							enhancer_str = enhancer_str + "NA,";
						}
					}
				}else{
							
					boolean[][] enhancer_spec_val = new boolean[dni_module.gene_region.length][dni_module.chip_tissue.length];
					
					for(int j = 0; j < enhancer_spec.length; j++){//region

						for(int l = 0; l < enhancer_spec[j].length; l++){//start
							
							for(int k = 0; k < enhancer_spec[j][l].length; k++){//tissue
							
								if(enhancer_spec[j][l][k].booleanValue()){
									enhancer_spec_val[j][k] = true;
									
								}

							}
						}
					}
					
					
					
					for(int j = 0; j < enhancer_spec_val.length; j++){//region

						for(int l = 0; l < enhancer_spec_val[j].length; l++){//tissue
						
							enhancer_str = enhancer_str + enhancer_spec_val[j][l] + ",";
						}
					}
				}//not null
				
				System.out.println("enhancer input:" + gene + "," + homologue_module[w][q] + "," + pcsd_str +  enhancer_str + homologue_path[w] );
					
				
				

				System.out.println("final input:" + gene + "," + homologue_module[w][q] + "," + pcsd_str + width+ tissue_str + enhancer_str + homologue_path[w] );
		  	}
	  		
		}// q lines
	
	
	
	
}//method

void test_init_ocr(){
String header = "";
	for(int z = 0; z < dni_module.type.length; z++){
		
		for(int x = 0; x < dni_module.gene_region.length; x++){

			for(int c = 0; c < dni_module.motif_list.length; c++){

				for(int v = 0; v < dni_module.strand_list.length; v++){

					header= header + dni_module.type[z] + "_" + dni_module.gene_region[x] + "_" + dni_module.motif_list[c] + "_" + dni_module.strand_list[v] + ",";

				}
			}
		}
	}

	

	//for(int i = 0; i < PO_GO.po_cluster.length; i++){

		//for(int e = 0; e < PO_GO.po_cluster_po_num[i].length; e++){

			//for(int r = 0; r < PO_GO.cluster_genes_po[i][e].size(); r++){

				//String gene = PO_GO.cluster_genes_po[i][e].get(r);
				String gene = "AT1G02230";
				LinkedList<String> temp_start = dni_module.peak_start.get(gene);
				
					//get gene motif, tf, LinkedHashMap<String,Boolean[]> motif_acc_to_domain_bool;LinkedHashMap<String,Boolean[]> domain_to_motif_acc_bool;, LinkedHashMap<String,LinkedList<String>> motif_acc_to_cluster LinkedHashMap<String,LinkedList<String>> domain_to_cluster;;
					
					Integer[][][][][] cur_motif = dni_module.gene_to_region_motif.get(gene);//motif_val[type_int][region_int][motif_int][strand_int]

					String write_motif = "";

					if(cur_motif != null){
					
						for(int z = 0; z < cur_motif.length; z++){
			
							for(int x = 0; x < cur_motif[z].length; x++){

								for(int c = 0; c < cur_motif[z][x].length; c++){

									for(int v = 0; v < cur_motif[z][x][c].length; v++){
										for(int w = 0; w < cur_motif[z][x][c][v].length; w++){

											write_motif = write_motif + cur_motif[z][x][c][v][w].intValue() + ",";
										}

									}
								}
							}
						}
					}else{

						for(int z = 0; z < dni_module.type.length; z++){
		
							for(int x = 0; x < dni_module.gene_region.length; x++){
for(int w = 0; w < temp_start.size(); w++){
								
									for(int c = 0; c < dni_module.motif_list.length; c++){

										for(int v = 0; v < dni_module.strand_list.length; v++){

											write_motif = write_motif + "null,";

										}
									}
								}
							}
						}
					}

					

					System.out.println("dni ocr:" + gene + "," + write_motif);
			//}
		//}
	//}//i
}

//test ddi
void test_ddi() throws IOException{

	for(int i = 0; i <PO_GO.cluster_genes_po.length; i++){
	
		for(int j =0; j < PO_GO.cluster_genes_po[i].length; j++){
		
			for(int c =0; c < PO_GO.cluster_genes_po[i][j].size(); c++){
		
				String gene = PO_GO.cluster_genes_po[i][j].get(c);
				LinkedList<String> pfams = ddi_interface_anno_module.atxg_to_pfam.get(gene);
				
				if(pfams != null){
					boolean[] sum_cat_stat_only_one = new boolean[ddi_interface_anno_module.cat_stat.length];
					boolean[] sum_cat_stat_only_multibp = new boolean[ddi_interface_anno_module.cat_stat.length];
					boolean[] sum_cat_stat_one_multibp = new boolean[ddi_interface_anno_module.cat_stat.length];
					
					//String[] aa_list = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};
					int[] max_aa = new int[ddi_interface_anno_module.aa_list.length];
					int[] min_aa = new int[ddi_interface_anno_module.aa_list.length];
					int sum_dimer_topo=0;
					
					if(pfams.size() < 5 && pfams.size() >=1){
					
						for(int k =0; k < pfams.size(); k++){
						
							String pfam1 = pfams.get(k);
							
							//System.out.println("pfam1:" + pfam1);
							
							if(pfam1 != null){
							
								
							
								String domain1 = ddi_interface_anno_module.pfam_to_domain.get(pfam1);
								//System.out.println("domain1:" + domain1);
								LinkedList<LinkedList<String>> all_sub_names = ddi_interface_anno_module.pfam_to_sub_folder.get(pfam1);
								
								if(all_sub_names != null){
								
									
								
									for(int l = 0; l < all_sub_names.size(); l++){
																		LinkedList<String> cur_folder_name = all_sub_names.get(l);	
										String query = "";
										for(int a = 0; a < cur_folder_name.size(); a++){
										
											query = query + cur_folder_name.get(a) + "@";						
										}
										
										query = query.substring(0,query.length()-1);
									//System.out.println("query:" + query);	//System.out.println("get_cat_stat_only_multibp(query):" + query);
										boolean[] cat_stat_only_multibp =  ddi_interface_anno_module.get_cat_stat_only_multibp(query);
										//System.out.println("return get_cat_stat_only_multibp(query):" + Arrays.toString(cat_stat_only_multibp));
										
										if(cat_stat_only_multibp == null){
										
										
										}else{
										
											for(int g =0; g < cat_stat_only_multibp.length; g++){
											
												if(cat_stat_only_multibp[g]){
												
													sum_cat_stat_only_multibp[g] = true;
												}
											}
										}
										
										//System.out.println("get_cat_stat_one_multibp(query):" + query);
										boolean[] cat_stat_one_multibp =  ddi_interface_anno_module.get_cat_stat_one_multibp(query);
										
										//System.out.println("return get_cat_stat_one_multibp(query):" + query + ":" + Arrays.toString(cat_stat_one_multibp));
										
										if(cat_stat_one_multibp == null){
										
										
										}else{
										
											for(int g =0; g < cat_stat_one_multibp.length; g++){
											
												if(cat_stat_one_multibp[g]){
												
													sum_cat_stat_one_multibp[g] = true;
												}
											}
										}
										
									}
									
									
									
									////System.out.println("get_cat_stat_only_onebp(pfam1):" + pfam1);
								}// all names not null
								
								boolean[] cat_stat_only_one =  ddi_interface_anno_module.get_cat_stat_only_onebp(pfam1);
								//System.out.println("return get_cat_stat_only_onebp(pfam1):" + Arrays.toString(cat_stat_only_one));
								
								if(cat_stat_only_one == null){
								
								
								}else{
								
									for(int l =0; l < cat_stat_only_one.length; l++){
									
										if(cat_stat_only_one[l]){
										
											sum_cat_stat_only_one[l] = true;
										}
									}
								}
								
								
								if(domain1 != null){
							
									for(int l =0; l < pfams.size(); l++){
									
										String partner_domain = ddi_interface_anno_module.pfam_to_domain.get(pfams.get(l));
										
										if(partner_domain != null){
											
											int dimer_topo_num = ddi_interface_anno_module.map_interface_topo_num(domain1 + "@" + partner_domain);
											sum_dimer_topo = sum_dimer_topo + dimer_topo_num;
											if(dimer_topo_num != 0){
											
												int[] sum_per_aa_topo_num = ddi_interface_anno_module.map_inter_topo(domain1 + "@" + partner_domain);
											
												if(sum_per_aa_topo_num != null){
													int max = sum_per_aa_topo_num[0];
													int max_ind = 0;
													
													int min = sum_per_aa_topo_num[0];
													int min_ind = 0;
													
													for(int z =1; z < sum_per_aa_topo_num.length; z++){
													
														if(max < sum_per_aa_topo_num[z]){
														
															max = sum_per_aa_topo_num[z];
															max_ind = z;
														}
														
														
														if(min> sum_per_aa_topo_num[z]){
														
															min = sum_per_aa_topo_num[z];
															min_ind = z;
														}
													}//for
													
													max_aa[max_ind] = max_aa[max_ind] + max;
													min_aa[min_ind] = min_aa[min_ind] + min;
													
												}//not null
													
													
												
											}//not null
											
										}//if not null
										
									}
								}
							}//if not null
						}//k
					}else{
					
						System.out.println("pfams >= 5:" + gene + ":" + pfams.size());
					
					}
					
					//write
					String cat_stat_write = gene + "," +PO_GO.po_cluster_po_num[i][j] + PO_GO.po_cluster_po_num[i][j] + ",";
					
					for(int z =0; z < sum_cat_stat_only_one.length; z++){
					
						cat_stat_write = cat_stat_write + sum_cat_stat_only_one[z] + ",";
					} 
					System.out.println("sum_cat_stat_only_one:" + cat_stat_write + sum_dimer_topo + "," + PO_GO.po_cluster[i]);
					
					//write
					cat_stat_write = gene + "," + PO_GO.po_cluster_po_num[i][j] + ",";
					
					for(int z =0; z < sum_cat_stat_only_multibp.length; z++){
					
						cat_stat_write = cat_stat_write + sum_cat_stat_only_multibp[z] + ",";
					} 
					System.out.println("sum_cat_stat_only_multibp:" + cat_stat_write + sum_dimer_topo + "," + PO_GO.po_cluster[i]);
					
					cat_stat_write = gene + "," + PO_GO.po_cluster_po_num[i][j] + ",";
					
					for(int z =0; z < sum_cat_stat_one_multibp.length; z++){
					
						cat_stat_write = cat_stat_write + sum_cat_stat_one_multibp[z] + ",";
					} 
					System.out.println("sum_cat_stat_one_multibp:" + cat_stat_write + sum_dimer_topo + "," + PO_GO.po_cluster[i]);
					
					//write max_aa
					String max_aa_write = gene + "," + PO_GO.po_cluster_po_num[i][j] + ",";
					
					for(int z =0; z < max_aa.length; z++){
					
						max_aa_write = max_aa_write +max_aa[z] + ",";
					} 
					System.out.println("max_aa_write:" + max_aa_write + sum_dimer_topo + "," + PO_GO.po_cluster[i]);
					
					String min_aa_write = gene + "," + PO_GO.po_cluster_po_num[i][j] + ",";
					
					for(int z =0; z < min_aa.length; z++){
					
						min_aa_write = min_aa_write + min_aa[z] + ",";
					} 
					System.out.println("min_aa_write:" + min_aa_write + sum_dimer_topo + "," + PO_GO.po_cluster[i]);
				}//if pfams not null
			
			}//size
			
		}//
	}//
				
	/*for(int i = 0; i <PO_GO.cluster_genes_slim.length; i++){
	
		for(int j =0; j < PO_GO.cluster_genes_slim[i].length; j++){
		
			for(int c =0; c < PO_GO.cluster_genes_slim[i][j].size(); c++){
		
				String gene = PO_GO.cluster_genes_slim[i][j].get(c);
				LinkedList<String> pfams = ddi_interface_anno_module.atxg_to_pfam.get(gene);
				
				if(pfams != null){
					boolean[] sum_cat_stat_only_one = new boolean[ddi_interface_anno_module.cat_stat.length];
					boolean[] sum_cat_stat_only_multibp = new boolean[ddi_interface_anno_module.cat_stat.length];
					boolean[] sum_cat_stat_one_multibp = new boolean[ddi_interface_anno_module.cat_stat.length];
					
					//String[] aa_list = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};
					int[] max_aa = new int[ddi_interface_anno_module.aa_list.length];
					int[] min_aa = new int[ddi_interface_anno_module.aa_list.length];
					int sum_dimer_toslim=0;
					
					if(pfams.size() < 5 && pfams.size() >=1){
					
						for(int k =0; k < pfams.size(); k++){
						
							String pfam1 = pfams.get(k);
							
							//System.out.println("pfam1:" + pfam1);
							
							if(pfam1 != null){
							
								
							
								String domain1 = ddi_interface_anno_module.pfam_to_domain.get(pfam1);
								//System.out.println("domain1:" + domain1);
								LinkedList<LinkedList<String>> all_sub_names = ddi_interface_anno_module.pfam_to_sub_folder.get(pfam1);
								
								if(all_sub_names != null){
								
									
								
									for(int l = 0; l < all_sub_names.size(); l++){
																		LinkedList<String> cur_folder_name = all_sub_names.get(l);	
										String query = "";
										for(int a = 0; a < cur_folder_name.size(); a++){
										
											query = query + cur_folder_name.get(a) + "@";						
										}
										
										query = query.substring(0,query.length()-1);
										//System.out.println("query:" + query);	//System.out.println("get_cat_stat_only_multibp(query):" + query);
										boolean[] cat_stat_only_multibp =  ddi_interface_anno_module.get_cat_stat_only_multibp(query);
										//System.out.println("return get_cat_stat_only_multibp(query):" + Arrays.toString(cat_stat_only_multibp));
										
										if(cat_stat_only_multibp == null){
										
										
										}else{
										
											for(int g =0; g < cat_stat_only_multibp.length; g++){
											
												if(cat_stat_only_multibp[g]){
												
													sum_cat_stat_only_multibp[g] = true;
												}
											}
										}
										
										//System.out.println("get_cat_stat_one_multibp(query):" + query);
										boolean[] cat_stat_one_multibp =  ddi_interface_anno_module.get_cat_stat_one_multibp(query);
										
										//System.out.println("return get_cat_stat_one_multibp(query):" + query + ":" + Arrays.toString(cat_stat_one_multibp));
										
										if(cat_stat_one_multibp == null){
										
										
										}else{
										
											for(int g =0; g < cat_stat_one_multibp.length; g++){
											
												if(cat_stat_one_multibp[g]){
												
													sum_cat_stat_one_multibp[g] = true;
												}
											}
										}
										
									}
									
									
									
									//System.out.println("get_cat_stat_only_onebp(pfam1):" + pfam1);
								}// all names not null
								
								boolean[] cat_stat_only_one =  ddi_interface_anno_module.get_cat_stat_only_onebp(pfam1);
								//System.out.println("return get_cat_stat_only_onebp(pfam1):" + Arrays.toString(cat_stat_only_one));
								
								if(cat_stat_only_one == null){
								
								
								}else{
								
									for(int l =0; l < cat_stat_only_one.length; l++){
									
										if(cat_stat_only_one[l]){
										
											sum_cat_stat_only_one[l] = true;
										}
									}
								}
								
								
								if(domain1 != null){
							
									for(int l =0; l < pfams.size(); l++){
									
										String partner_domain = ddi_interface_anno_module.pfam_to_domain.get(pfams.get(l));
										
										if(partner_domain != null){
											
											int dimer_toslim_num = ddi_interface_anno_module.map_interface_topo_num(domain1 + "@" + partner_domain);
											sum_dimer_toslim = sum_dimer_toslim + dimer_toslim_num;
											if(dimer_toslim_num != 0){
											
												int[] sum_per_aa_toslim_num = ddi_interface_anno_module.map_inter_topo(domain1 + "@" + partner_domain);
											
												if(sum_per_aa_toslim_num != null){
													int max = sum_per_aa_toslim_num[0];
													int max_ind = 0;
													
													int min = sum_per_aa_toslim_num[0];
													int min_ind = 0;
													
													for(int z =1; z < sum_per_aa_toslim_num.length; z++){
													
														if(max < sum_per_aa_toslim_num[z]){
														
															max = sum_per_aa_toslim_num[z];
															max_ind = z;
														}
														
														
														if(min> sum_per_aa_toslim_num[z]){
														
															min = sum_per_aa_toslim_num[z];
															min_ind = z;
														}
													}//for
													
													max_aa[max_ind] = max_aa[max_ind] + max;
													min_aa[min_ind] = min_aa[min_ind] + min;
													
												}//not null
													
													
												
											}//not null
											
										}//if not null
										
									}
								}
							}//if not null
						}//k
					}else{
					
						//System.out.println("pfams >= 5:" + gene + ":" + pfams.size());
					
					}
					
					//write
					String cat_stat_write = gene + "," +PO_GO.slim_cluster_slim_num[i][j] + PO_GO.slim_cluster_slim_num[i][j] + ",";
					
					for(int z =0; z < sum_cat_stat_only_one.length; z++){
					
						cat_stat_write = cat_stat_write + sum_cat_stat_only_one[z] + ",";
					} 
					System.out.println("sum_cat_stat_only_one:" + cat_stat_write + sum_dimer_toslim + "," + PO_GO.slim_cluster[i]);
					
					//write
					cat_stat_write = gene + "," + PO_GO.slim_cluster_slim_num[i][j] + ",";
					
					for(int z =0; z < sum_cat_stat_only_multibp.length; z++){
					
						cat_stat_write = cat_stat_write + sum_cat_stat_only_multibp[z] + ",";
					} 
					System.out.println("sum_cat_stat_only_multibp:" + cat_stat_write + sum_dimer_toslim + "," + PO_GO.slim_cluster[i]);
					
					cat_stat_write = gene + "," + PO_GO.slim_cluster_slim_num[i][j] + ",";
					
					for(int z =0; z < sum_cat_stat_one_multibp.length; z++){
					
						cat_stat_write = cat_stat_write + sum_cat_stat_one_multibp[z] + ",";
					} 
					System.out.println("sum_cat_stat_one_multibp:" + cat_stat_write + sum_dimer_toslim + "," + PO_GO.slim_cluster[i]);
					
					//write max_aa
					String max_aa_write = gene + "," + PO_GO.slim_cluster_slim_num[i][j] + ",";
					
					for(int z =0; z < max_aa.length; z++){
					
						max_aa_write = max_aa_write +max_aa[z] + ",";
					} 
					System.out.println("max_aa_write:" + max_aa_write + sum_dimer_toslim + "," + PO_GO.slim_cluster[i]);
					
					String min_aa_write = gene + "," + PO_GO.slim_cluster_slim_num[i][j] + ",";
					
					for(int z =0; z < min_aa.length; z++){
					
						min_aa_write = min_aa_write + min_aa[z] + ",";
					} 
					System.out.println("min_aa_write:" + min_aa_write + sum_dimer_toslim + "," + PO_GO.slim_cluster[i]);
				}//if pfams not null
			
			}//size
			
		}//
	}//
	*/			
					
}

//////test demonstration functions comes here///

void test_compound() throws IOException{

	//System.out.println("map_pheno_to_path():" );
	
	for(int w = 0; w <PO_GO.cluster_genes_po.length; w++){
	
		for(int q =0; q < PO_GO.cluster_genes_po[w].length; q++){
		
			for(int c =0; c < PO_GO.cluster_genes_po[w][q].size(); c++){
		
				String gene = PO_GO.cluster_genes_po[w][q].get(c);
				LinkedList<String> cur_prod_enz = pathway_anal_module.pmn_acc_to_product[0].get(gene);
				
				if(cur_prod_enz != null){
				
					for(int g = 0; g < cur_prod_enz.size(); g++){
				
						String enzrxn_id = pathway_anal_module.product_to_enzrxn_id[0].get(cur_prod_enz.get(g));
						//System.out.println("183 enzrxn_id:" + gene + ":" +enzrxn_id);
		 				Boolean[] cur_rxn_type =  pathway_anal_module.enzrxn_id_to_rxn_type[0].get(enzrxn_id);

						String write_rxn_type = "";
						if(cur_rxn_type != null){
							
							for(int z = 0; z < cur_rxn_type.length; z++){

								write_rxn_type = write_rxn_type + cur_rxn_type[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.rxn_type.length; z++){
								write_rxn_type = write_rxn_type  + "null,";
				
							}

						}
						String cur_pwy = pathway_anal_module.enzrxn_id_to_pwy[0].get(enzrxn_id);

						String write_pwy = "";
						if(cur_pwy != null){
							
							write_pwy =cur_pwy + ",";
						}else{
							write_pwy = "null,";
				
							

						}

						Boolean[] cur_substrate_form_ele = pathway_anal_module.enzrxn_id_to_substrate_form_ele[0].get(enzrxn_id);

						String write_substrate_form_ele = "";
						if(cur_substrate_form_ele != null){
							
							for(int z = 0; z < cur_substrate_form_ele.length; z++){

								write_substrate_form_ele = write_substrate_form_ele + cur_substrate_form_ele[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.formula_ele.length; z++){
								write_substrate_form_ele = write_substrate_form_ele  + "null,";
				
							}

						}

						Boolean[] cur_substrate_struc_group = pathway_anal_module.enzrxn_id_to_substrate_struc_group[0].get(enzrxn_id);

						String write_substrate_struc_group = "";
						if(cur_substrate_struc_group != null){
							
							for(int z = 0; z < cur_substrate_struc_group.length; z++){

								write_substrate_struc_group = write_substrate_struc_group + cur_substrate_struc_group[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.struc_group.length; z++){
								write_substrate_struc_group = write_substrate_struc_group  + "null,";
				
							}

						}

						Boolean[] cur_substrate_type_group = pathway_anal_module.enzrxn_id_to_substrate_type_group[0].get(enzrxn_id);

						String write_substrate_type_group = "";
						if(cur_substrate_type_group != null){
							
							for(int z = 0; z < cur_substrate_type_group.length; z++){

								write_substrate_type_group = write_substrate_type_group + cur_substrate_type_group[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.comp_type_group.length; z++){
								write_substrate_type_group = write_substrate_type_group  + "null,";
				
							}

						}

						Boolean[] cur_product_form_ele = pathway_anal_module.enzrxn_id_to_product_form_ele[0].get(enzrxn_id);
						String write_product_form_ele = "";
						if(cur_product_form_ele != null){
							
							for(int z = 0; z < cur_product_form_ele.length; z++){

								write_product_form_ele = write_product_form_ele + cur_product_form_ele[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.formula_ele.length; z++){
								write_product_form_ele = write_product_form_ele  + "null,";
				
							}

						}

						Boolean[] cur_product_struc_group = pathway_anal_module.enzrxn_id_to_product_struc_group[0].get(enzrxn_id);

						String write_product_struc_group = "";
						if(cur_product_struc_group != null){
							
							for(int z = 0; z < cur_product_struc_group.length; z++){

								write_product_struc_group = write_product_struc_group + cur_product_struc_group[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.struc_group.length; z++){
								write_product_struc_group = write_product_struc_group  + "null,";
				
							}

						}

						
						Boolean[] cur_product_type_group = pathway_anal_module.enzrxn_id_to_product_type_group[0].get(enzrxn_id);

						String write_product_type_group = "";
						if(cur_product_type_group != null){
							
							for(int z = 0; z < cur_product_type_group.length; z++){

								write_product_type_group = write_product_type_group + cur_product_type_group[z].booleanValue() + ",";
							}
						}else{
							for(int z = 0; z < pathway_anal_module.comp_type_group.length; z++){
								write_product_type_group = write_product_type_group  + "null,";
				
							}

						}


						LinkedList<String> left_compounds = pathway_anal_module.enzrxn_id_to_substrate_compound[0].get(enzrxn_id);
						
						
						double left_avr = 0.0;
						if(left_compounds == null){
							//System.out.println("error left compound id:" +"null");

						}else{
							if(left_compounds.size() == 1){
				
								if(left_compounds.get(0) != null){

										if(pathway_anal_module.compound_id_to_mw[0].get(left_compounds.get(0)) != null){
										
										//System.out.println("error compound id:" + left_compounds.get(0));	
										left_avr = pathway_anal_module.compound_id_to_mw[0].get(left_compounds.get(0)).doubleValue();	
									}		
								}
									

							}else{

								double[] val = new double[left_compounds.size()];

								for(int f = 0; f < left_compounds.size(); f++){

									if(left_compounds.get(f) != null){
										//System.out.println("error compound id 2:" + left_compounds.get(f));
										if(pathway_anal_module.compound_id_to_mw[0].get(left_compounds.get(f)) != null){
											val[f] = pathway_anal_module.compound_id_to_mw[0].get(left_compounds.get(f)).doubleValue();
										}
									}

								}
								DescriptiveStatistics desc = new DescriptiveStatistics(val);
								left_avr=new Double(desc.getMean());
								


							}
						}

						LinkedList<String> right_compounds = pathway_anal_module.enzrxn_id_to_substrate_compound[0].get(enzrxn_id);

						double right_avr = 0.0;

						if(right_compounds == null){

							//System.out.println("error right compound id:" +"null");
						}else{

							if(right_compounds.size() == 1){

								if(right_compounds.get(0) != null){

									if(pathway_anal_module.compound_id_to_mw[0].get(right_compounds.get(0)) != null){
										//System.out.println("error compound id 3:" + right_compounds.get(0));
										right_avr = pathway_anal_module.compound_id_to_mw[0].get(right_compounds.get(0)).doubleValue();
									}
								}

							}else{

								double[] val = new double[right_compounds.size()];

								for(int f = 0; f < right_compounds.size(); f++){

									if(right_compounds.get(f) != null){

									if(pathway_anal_module.compound_id_to_mw[0].get(right_compounds.get(f)) != null){
										//System.out.println("error compound id 3:" + right_compounds.get(f));
										val[f] = pathway_anal_module.compound_id_to_mw[0].get(right_compounds.get(f)).doubleValue();
										}
									}

								}
								DescriptiveStatistics desc = new DescriptiveStatistics(val);
								right_avr=new Double(desc.getMean());
								


							}
						}

						
						if(write_pwy.equals(",")){
						
							write_pwy = "NA,";
						}
							
						
						System.out.println("datasheet large pmn_test po_go:" + enzrxn_id + "," + cur_prod_enz.get(g) + "," + gene + "," +PO_GO.po_cluster_po_num[w][q] + "," +write_rxn_type + write_pwy + write_substrate_form_ele + write_substrate_struc_group + write_substrate_type_group + write_product_form_ele + write_product_struc_group + write_product_type_group + left_avr + "," + right_avr +"," + PO_GO.po_cluster[w]);
					}//for g
				}//not null
				
			}//q
		}//w
	}
}//method



////////////////////////////////
static class dni_module{




static LinkedHashMap<String, Integer[]> peak_size;
static LinkedHashMap<String, LinkedList<Range<Integer>>> peak_range;

static LinkedHashMap<String, LinkedList<String>> peak_start;
//get only high specificity regions. exclude low
void get_how_many_peak() throws IOException{

	Pattern space_pattern = Pattern.compile("\\s+");

	peak_size = new LinkedHashMap();
	peak_start = new LinkedHashMap();
	peak_range = new LinkedHashMap();

	
	String path = "data/arabidopsis_thaliana.final.annotatePeak.peak.start";
/*1	1190	1993	.	.	.	0.0281692132085616	AT1G01010
1	10877	11591	.	.	.	0.0296482466761526	AT1G01020
1	16423	16936	.	.	.	0.0284302785684615	AT1G01030
1	19366	19583	.	.	.	0.030104200112554	AT1G01040
1	39748	40640	.	.	.	0.0287764769777116	AT1G01060

*/
	List<String> lines  = FileUtils.readLines(new File(path));

        LinkedList<String> input = new LinkedList();

	Integer[] temp_size = new Integer[3];//freq in each of <500, >500<2000,>2000

	for(int i = 0; i < 3; i++){

		temp_size[i] = new Integer(0);
	}
        LinkedList<String> temp_pos = new LinkedList();
	String id = "";
	String pos = "";
	String pos1 = "";
	String pos2 = "";
	Integer size = new Integer(0);

	for(int q= 1; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
		String[] pos_split = space_pattern.split(one_line);
		String next_line = "";
		String[] next_split = new String[8];
		if(q != lines.size()-1){
			next_line = lines.get(q+1).trim();
			next_split = space_pattern.split(next_line);
		}

		if(!pos_split[7].equals(next_split[7]) || q == lines.size()-1){

			System.out.println("1165:"+one_line);

			if(q != lines.size()-1){
				//put into map
				System.out.println("1169:"+id + ":" + Arrays.toString(temp_size) + ":" + Arrays.toString(temp_pos.toArray()));
				//peak_size: 
				peak_size.put(id,temp_size);
				//start pos and range
				LinkedList<String> temp_pos_start = new LinkedList();
				LinkedList<Range<Integer>> temp_pos_range = new LinkedList();
				for(int f = 0; f < temp_pos.size(); f++){
				
					String position = temp_pos.get(f).trim();
					String position1 = position.substring(0, position.indexOf("_"));
					String position2 = position.substring(position.indexOf("_")+1,position.length());
					
					if(!temp_pos_start.contains(position1)){
					
						temp_pos_start.add(position1);
					}
					
					
					Integer start = new Integer(position1);
					Integer end = new Integer(position2);	
					Range<Integer> range = Range.between(start,end);
					temp_pos_range.add(range);		
					
				}
				peak_start.put(id,temp_pos_start);
				peak_range.put(id,temp_pos_range);
			
		

				//reset list
				temp_pos = new LinkedList();
				temp_size = new Integer[3];
				for(int i = 0; i < 3; i++){

					temp_size[i] = new Integer(0);
				}

			}else{
				pos1 = pos_split[1];
				pos2 = pos_split[2];
				
				if(!temp_pos.contains(pos1+ "_" + pos2)){
					
					if(temp_pos != null){
						temp_pos.add(pos1+ "_" + pos2);
					}
				}
				System.out.println("858:"+id+ ":" + Arrays.toString(temp_pos.toArray()));
				size = new Integer(pos2) - new Integer(pos1); 

				if(size >1 && size < 500){

					temp_size[0] = temp_size[0] +1;
				}else if(size >= 500 && size < 2000){

					temp_size[1] = temp_size[1] +1;
				}else if(size >2000){

					temp_size[2] = temp_size[2] +1;
				}
				//peak_size: 
				peak_size.put(id,temp_size);
				
				//start position and range
				LinkedList<String> temp_pos_start = new LinkedList();
				LinkedList<Range<Integer>> temp_pos_range = new LinkedList();
				for(int f = 0; f < temp_pos.size(); f++){
				
					String position = temp_pos.get(f).trim();
					String position1 = position.substring(0, position.indexOf("_"));
					String position2 = position.substring(position.indexOf("_")+1,position.length());
					
					if(!temp_pos_start.contains(position1)){
					
						temp_pos_start.add(position1);
					}
					
					
					Integer start = new Integer(position1);
					Integer end = new Integer(position2);	
					Range<Integer> range = Range.between(start,end);
					temp_pos_range.add(range);		
					
				}
				peak_start.put(id,temp_pos_start);
				peak_range.put(id,temp_pos_range);
			
				System.out.println("854:"+temp_pos.size() + ":" + Arrays.toString(temp_pos.toArray()));
			}

			id = pos_split[7];
			pos1 = pos_split[1];
			pos2 = pos_split[2];
			
			if(!temp_pos.contains(pos1+ "_" + pos2)){
				
				if(temp_pos != null){
					temp_pos.add(pos1+ "_" + pos2);
				}
			}
			System.out.println("861:"+id+ ":" + Arrays.toString(temp_pos.toArray()));
			size = new Integer(pos2) - new Integer(pos1); 

			if(size >1 && size < 500){

				temp_size[0] = temp_size[0] +1;
			}else if(size >= 500 && size < 2000){

				temp_size[1] = temp_size[1] +1;
			}else if(size >2000){

				temp_size[2] = temp_size[2] +1;
			}

			

		}else{

			id = pos_split[7];
			pos1 = pos_split[1];
			pos2 = pos_split[2];
			
			if(!temp_pos.contains(pos1+ "_" + pos2)){
				if(temp_pos != null){
				
					temp_pos.add(pos1+ "_" + pos2);
				}
			}
			size = new Integer(pos2) - new Integer(pos1); 

			if(size >1 && size < 500){

				temp_size[0] = temp_size[0] +1;
			}else if(size >= 500 && size < 2000){

				temp_size[1] = temp_size[1] +1;
			}else if(size >2000){

				temp_size[2] = temp_size[2] +1;
			}

		}//else

		
	}//q
	
}//method

static LinkedHashMap<String, Integer[]> peak_size_spec;
static LinkedHashMap<String, LinkedList<String>> peak_start_spec;
//get only high specificity regions. exclude low
void get_how_many_peak_spec() throws IOException{

	Pattern space_pattern = Pattern.compile("\\s+");

	peak_size_spec = new LinkedHashMap();
	peak_start_spec = new LinkedHashMap();

	
	String path = "data/specificity.all.peak_start.sort.uniq2";
	
/*1	1250	1843	.	.	.	.	AT1G01010
1	1957	3961	.	.	.	.	AT1G01010
1	8087	9108	.	.	.	.	AT1G01020


*/
	List<String> lines  = FileUtils.readLines(new File(path));

        LinkedList<String> input = new LinkedList();

	Integer[] temp_size = new Integer[3];//freq in each of <500, >500<2000,>2000

	for(int i = 0; i < 3; i++){

		temp_size[i] = new Integer(0);
	}
        LinkedList<String> temp_pos = new LinkedList();
	String id = "";
	String pos = "";
	String pos1 = "";
	String pos2 = "";
	Integer size = new Integer(0);

	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
		String[] pos_split = space_pattern.split(one_line);
		String next_line = "";
		String[] next_split = new String[8];
		if(q != lines.size()-1){
			next_line = lines.get(q+1).trim();
			next_split = space_pattern.split(next_line);
		}

		if(!pos_split[7].equals(next_split[7]) || q == lines.size()-1){

			System.out.println("1165:"+one_line);

			if(q != lines.size()-1){
				//put into map
				System.out.println("1169:"+id + ":" + Arrays.toString(temp_size) + ":" + Arrays.toString(temp_pos.toArray()));
				//peak_size_spec: 
				peak_size_spec.put(id,temp_size);
				//start pos
				peak_start_spec.put(id,temp_pos);

				//reset list
				temp_pos = new LinkedList();
				temp_size = new Integer[3];
				for(int i = 0; i < 3; i++){

					temp_size[i] = new Integer(0);
				}

			}else{
				pos1 = pos_split[1];
				pos2 = pos_split[2];
				
				if(!temp_pos.contains(pos1)){
					
					if(temp_pos != null){
						temp_pos.add(pos1);
					}
				}
				System.out.println("858:"+id+ ":" + Arrays.toString(temp_pos.toArray()));
				size = new Integer(pos2) - new Integer(pos1); 

				if(size >1 && size < 500){

					temp_size[0] = temp_size[0] +1;
				}else if(size >= 500 && size < 2000){

					temp_size[1] = temp_size[1] +1;
				}else if(size >2000){

					temp_size[2] = temp_size[2] +1;
				}
				//peak_size_spec: 
				peak_size_spec.put(id,temp_size);
				//start pos
				peak_start_spec.put(id,temp_pos);
				System.out.println("854:"+temp_pos.size() + ":" + Arrays.toString(temp_pos.toArray()));
			}

			id = pos_split[7];
			pos1 = pos_split[1];
			pos2 = pos_split[2];
			
			if(!temp_pos.contains(pos1)){
				
				if(temp_pos != null){
					temp_pos.add(pos1);
				}
			}
			System.out.println("861:"+id+ ":" + Arrays.toString(temp_pos.toArray()));
			size = new Integer(pos2) - new Integer(pos1); 

			if(size >1 && size < 500){

				temp_size[0] = temp_size[0] +1;
			}else if(size >= 500 && size < 2000){

				temp_size[1] = temp_size[1] +1;
			}else if(size >2000){

				temp_size[2] = temp_size[2] +1;
			}

			

		}else{

			id = pos_split[7];
			pos1 = pos_split[1];
			pos2 = pos_split[2];
			
			if(!temp_pos.contains(pos1)){
				if(temp_pos != null){
				
					temp_pos.add(pos1);
				}
			}
			size = new Integer(pos2) - new Integer(pos1); 

			if(size >1 && size < 500){

				temp_size[0] = temp_size[0] +1;
			}else if(size >= 500 && size < 2000){

				temp_size[1] = temp_size[1] +1;
			}else if(size >2000){

				temp_size[2] = temp_size[2] +1;
			}

		}//else

		
	}//q
	
}//method


static String[] type = {"Promoter","Enhancer"};
static String[] gene_region={"3’UTR","5’UTR","Exon","Intergenic","Intron"};

static String[] tf_list ={
"ABF1","ABF2","ABF3","ABF4","ABI3","abi4","ABI5","ABR1","AG","AGL1","AGL13","AGL15","AGL16","AGL27","AGL3","AGL42","AGL55","AGL6","AGL63","AHL12","AHL20","AHL25","AIB","AIL6","AIL7","ANL2","ANT","AP1","AP3","ARALYDRAFT_484486","ARALYDRAFT_493022","ARALYDRAFT_495258","ARALYDRAFT_496250","ARALYDRAFT_897773","ARF1","ARF10","ARF13","ARF14","ARF16","ARF18","ARF2","ARF25","ARF27","ARF29","ARF3","ARF34","ARF35","ARF36","ARF39","ARF4","ARF5","ARF7","ARF8","ARR1","ARR10","ARR101","ARR11","ARR14","ARR18","ARR2","ASIL2","ASR1","ASR3","AT1G11260","AT1G14240","AT1G14600","AT1G19000","AT1G19040","AT1G21920","AT1G29920","AT1G72740","AT1G74840","AT1G76870","AT2G38090","AT2G38300","AT2G40260","AT3G10030","AT3G10580","At3g11280","AT3G46070","AT5G04390","AT5G04760","AT5G05090","AT5G05550","AT5G05790","AT5G47660","AT5G56840","AT5G61620","ATHB-12","ATHB-13","ATHB-15","ATHB-16","ATHB-20","ATHB-23","ATHB-4","ATHB-40","ATHB-5","ATHB-51","ATHB-53","ATHB-6","ATHB-7","ATHB-9","ATHB-X","ATMYB31","BAM8","BEE2","BEH2","BEH3","BEH4","BHLH104","BHLH112","BHLH122","BHLH13","bHLH130","bHLH145","bHLH18","BHLH3","BHLH34","BHLH49","BHLH72","BHLH74","bHLH77","BHLH78","bHLH80","BIM1","BIM2","BIM3","BPC1","BPC5","BPC6","BPE","BRN2","BZIP11","BZIP16","BZIP18","BZIP2","BZIP28","BZIP3","BZIP30","BZIP42","BZIP43","BZIP44","BZIP48","BZIP52","BZIP53","BZIP60","BZIP63","BZIP68","BZIP69","bZIP910","bZIP911","BZR1","BZR2","CAMTA1","CAMTA2","CAMTA3","CCA1","CDC5","CDF2","CDF3","CDF5","CRF2","CRF4","DEAR3","DF1","DIV1","DOF1.5","DOF1.6","DOF1.7","DOF1.8","Dof2","DOF2.2","DOF2.4","DOF2.5","Dof3","DOF3.2","DOF3.4","DOF3.5","DOF3.6","DOF4.2","DOF4.3","DOF4.5","DOF5.1","DOF5.3","DOF5.4","DOF5.6","DOF5.7","DOF5.8","DPBF3","DREB1A","DREB1B","DREB1C","DREB1D","DREB1E","DREB1F","DREB1G","DREB2A","DREB2C","DREB2D","DREB2E","DREB2F","DREB2G","DYT1","E2FA","E2FC","E2FD","E2FE","EFM","EIL3","EIL4","EmBP-1","Enhancer","ERF003","ERF008","ERF010","ERF011","ERF012","ERF013","ERF014","ERF015","ERF017","ERF018","ERF019","ERF021","ERF023","ERF025","ERF027","ERF034","ERF035","ERF036","ERF037","ERF038","ERF039","ERF043","ERF054","ERF055","ERF057","ERF069","ERF073","ERF086","ERF087","ERF091","ERF094","ERF095","ERF096","ERF10","ERF104","ERF105","ERF109","ERF11","ERF112","ERF115","ERF118","ERF122","ERF13","ERF15","ERF1B","ERF2","ERF3","ERF4","ERF5","ERF6","ERF7","ERF8","ERF9","Exon","FaEOBII","FAR1","FHY3","FLC","FUS3","GAF1","Gam1","GATA10","GATA11","GATA12","GATA14","GATA15","GATA19","GATA20","GATA4","GATA6","GATA8","GATA9","GBF2","GBF3","GLYMA-06G314400","GLYMA-07G038400","GLYMA-08G357600","GLYMA-13G317000","Glyma19g26560.1","GRF4","GRF6","GRF9","GT-1","GT-2","GT-3a","GT-4","GTL1","HAT1","HAT2","HAT22","HAT5","HBI1","HDG1","HDG11","HDG7","HHO2","HHO3","HHO5","HHO6","HRS1","HSFA1B","HSFA1E","HSFA4A","HSFA6A","HSFA6B","HSFB2A","HSFB2B","HSFB3","HSFB4","HSFC1","HY5","HYH","IDD4","IDD5","IDD6","IDD7","Intergenic","JKD","JUB1","KAN1","KAN2","KAN4","KUA1","LBD13","LBD18","LEC2","LEP","LFY","LHY","LOB","LRL2","MA1186.1","MA1362.1","MA1682.1","MGP","MNB1A","MYB1","MYB10","MYB101","MYB105","MYB107","MYB108","MYB111","MYB113","MYB116","MYB118","MYB119","MYB121","MYB124","MYB13","MYB15","MYB17","MYB23","MYB24","MYB27","MYB3","MYB30","MYB33","MYB39","MYB3R1","MYB3R4","MYB3R5","MYB4","MYB40","MYB41","MYB43","MYB44","MYB46","MYB49","MYB51","MYB52","MYB55","MYB56","MYB57","MYB58","MYB59","MYB60","MYB61","MYB62","MYB63","MYB65","MYB67","MYB70","MYB73","MYB74","MYB77","MYB80","MYB81","MYB83","MYB88","MYB92","MYB93","MYB94","MYB96","MYB98","MYB99","myb.Ph3","MYC2","MYC3","MYC4","MYR2","NAC002","NAC004","NAC005","NAC007","NAC010","NAC011","NAC013","NAC016","NAC017","NAC018","NAC019","NAC020","NAC025","NAC028","NAC029","NAC031","NAC035","NAC037","NAC038","NAC043","NAC045","NAC046","NAC047","NAC050","NAC053","NAC054","NAC055","NAC057","NAC058","NAC062","NAC071","NAC073","NAC075","NAC076","NAC0761","NAC078","NAC079","NAC083","NAC087","NAC096","NAC098","NAC101","NAC103","NAC105","NAC68","NAC69","NAC92","NID1","NLP7","NTL8","NTL9","NUC","O11","O2","OJ1058_F05.8","OJ1581_H09.2","Os05g0497200","OsI_08196","OsRR22","P0510F09.23","PBF","PEND","PHL1","PHL11","PHL12","PHL2","PHL4","PHL7","PHYPADRAFT_140773","PHYPADRAFT_143875","PHYPADRAFT_153324","PHYPADRAFT_173530","PHYPADRAFT_182268","PHYPADRAFT_28324","PHYPADRAFT_38837","PHYPADRAFT_48267","PHYPADRAFT_64121","PHYPADRAFT_72483","PI","PIF1","PIF3","PIF4","PIF5","PLT1","POPTR_0002s00440g","RAMOSA1","RAP1","RAP2-1","RAP2-10","RAP2-11","RAP2-12","RAP2-3","RAP2-4","RAP2-6","RAP2-7","RAP2-9","RAV1","RAV2","RAX3","REF6","RIN","RVE1","RVE4","RVE5","RVE6","RVE7","RVE7L","RVE8","SEP1","SEP3","SGR5","SIZF2","SMB","SMZ","SOC1","SOL1","SPL1","SPL10","SPL11","SPL12","SPL13A","SPL14","SPL15","SPL3","SPL4","SPL5","SPL7","SPL8","SPL9","SPT","squamosa","SRM1","StBRC1","SVP","TB1","TCP1","TCP13","TCP14","TCP15","TCP16","TCP17","TCP19","TCP2","TCP20","TCP21","TCP22","TCP23","TCP24","TCP3","TCP4","TCP5","TCP7","TCP8","TCP9","TCX2","TCX3","TCX6","TEM1","TFLG2-Zm00001d042777","TGA1","TGA10","TGA1A","TGA2","TGA3","TGA4","TGA5","TGA6","TGA7","TGA9","TINY","TRB1","TRB2","TREE1","TRP1","TRP2","TRP5","TSAR1","TSAR2","TSO1","UNE10","VIP1","WIN1","WRKY1","WRKY11","WRKY12","WRKY14","WRKY15","WRKY17","WRKY18","WRKY2","WRKY20","WRKY21","WRKY22","WRKY23","WRKY24","WRKY25","WRKY26","WRKY27","WRKY28","WRKY29","WRKY3","WRKY30","WRKY31","WRKY33","WRKY38","WRKY40","WRKY42","WRKY43","WRKY45","WRKY46","WRKY47","WRKY48","WRKY50","WRKY53","WRKY55","WRKY57","WRKY59","WRKY6","WRKY60","WRKY62","WRKY63","WRKY65","WRKY7","WRKY70","WRKY71","WRKY75","WRKY8","YAB4","ZAT10","ZAT6","ZHD1","ZHD10","ZHD3","ZHD5","ZHD6","ZHD9","Zm00001d002364","Zm00001d005692","Zm00001d005892","Zm00001d015407","Zm00001d018571","Zm00001d020267","Zm00001d020595","Zm00001d024324","Zm00001d027846","Zm00001d031796","Zm00001d034298","Zm00001d035604","Zm00001d038683","Zm00001d044409","Zm00001d044785","Zm00001d049364","Zm00001d052229"
}; 

static String[] motif_list ={
"MA0001.2","MA0005.2","MA0008.3","MA0020.1","MA0021.1","MA0034.1","MA0053.1","MA0054.1","MA0064.1","MA0082.1","MA0096.1","MA0097.1","MA0110.3","MA0121.1","MA0123.1","MA0127.1","MA0128.1","MA0129.1","MA0548.2","MA0549.1","MA0550.2","MA0551.1","MA0552.1","MA0553.1","MA0554.1","MA0555.1","MA0556.1","MA0557.1","MA0558.1","MA0559.1","MA0560.1","MA0561.1","MA0562.1","MA0563.1","MA0564.1","MA0565.2","MA0566.1","MA0567.1","MA0568.1","MA0569.1","MA0570.2","MA0571.1","MA0573.1","MA0574.1","MA0575.1","MA0576.1","MA0577.3","MA0578.1","MA0579.1","MA0580.1","MA0581.1","MA0582.1","MA0583.1","MA0584.1","MA0585.1","MA0586.2","MA0587.1","MA0588.1","MA0589.1","MA0590.1","MA091","MA0930.2","MA0931.1","MA0932.1","MA0933.1","MA0934.1","MA0935.1","MA0936.1","MA0937.1","MA0938.2","MA0939.1","MA0940.1","MA0941.1","MA0942.1","MA0943.1","MA0944.1","MA0945.1","MA0946.1","MA0947.1","MA0948.1","MA0949.1","MA0950.1","MA0951.1","MA0952.1","MA0953.1","MA0954.2","MA0955.1","MA0956.1","MA0957.1","MA0958.1","MA0959.1","MA0960.1","MA0961.1","MA0962.1","MA0963.1","MA0964.2","MA0965.2","MA0966.1","MA0967.1","MA0968.2","MA0969.1","MA0970.1","MA0971.1","MA0972.1","MA0973.1","MA0974.2","MA0975.1","MA0976.2","MA0977.1","MA0978.1","MA0979.1","MA0980.2","MA0981.1","MA0982.1","MA0983.1","MA0984.1","MA0986.1","MA0987.1","MA0988.1","MA0989.1","MA0990.1","MA0992.2","MA0993.1","MA0994.2","MA0995.2","MA0997.1","MA0998.1","MA1000.2","MA1001.3","MA1004.1","MA1005.2","MA1006.1","MA1007.1","MA1008.1","MA1009.1","MA1010.1","MA1011.1","MA1012.1","MA1013.1","MA1014.1","MA1015.1","MA1016.1","MA1017.1","MA1018.1","MA1019.1","MA1020.1","MA1021.1","MA1022.1","MA1023.1","MA1024.1","MA1025.1","MA1026.2","MA1027.1","MA1028.1","MA1030.1","MA1031.1","MA1032.1","MA1033.1","MA1034.1","MA1035.1","MA1036.1","MA1037.1","MA1038.1","MA1039.1","MA1040.1","MA1041.1","MA1042.1","MA1043.1","MA1044.1","MA1045.1","MA1046.1","MA1047.2","MA1048.1","MA1049.1","MA1050.1","MA1051.1","MA1053.1","MA1054.1","MA1055.2","MA1056.1","MA1057.1","MA1058.1","MA1059.2","MA1060.1","MA1061.1","MA1062.2","MA1063.1","MA1064.1","MA1065.2","MA1066.1","MA1067.1","MA1068.2","MA1069.2","MA1070.2","MA1071.1","MA1074.1","MA1075.1","MA1076.2","MA1077.1","MA1078.1","MA1079.2","MA1080.1","MA1081.2","MA1083.2","MA1084.1","MA1085.2","MA1086.1","MA1087.2","MA1088.1","MA1089.1","MA1090.1","MA1091.1","MA1092.1","MA1093.1","MA1094.2","MA1095.1","MA1096.1","MA1097.1","MA1098.1","MA1156.1","MA1157.1","MA1158.1","MA1159.1","MA1160.1","MA1161.1","MA1162.1","MA1163.2","MA1164.1","MA1165.2","MA1166.1","MA1167.1","MA1168.1","MA1169.1","MA1170.1","MA1171.1","MA1172.1","MA1173.1","MA1174.1","MA1175.1","MA1176.1","MA1177.1","MA1178.2","MA1179.1","MA1180.1","MA1181.1","MA1182.1","MA1183.1","MA1184.1","MA1185.1","MA1186.1","MA1187.1","MA1188.1","MA1189.1","MA1190.1","MA1191.1","MA1192.1","MA1193.1","MA1194.1","MA1195.1","MA1196.1","MA1197.1","MA1198.1","MA1199.1","MA1201.1","MA1202.1","MA1203.1","MA1204.1","MA1205.1","MA1206.1","MA1207.1","MA1208.1","MA1209.1","MA1210.2","MA1211.1","MA1212.1","MA1213.2","MA1214.1","MA1215.1","MA1218.1","MA1219.2","MA1220.1","MA1221.1","MA1222.1","MA1223.1","MA1225.1","MA1226.1","MA1227.2","MA1228.1","MA1229.1","MA1230.1","MA1231.2","MA1232.1","MA1233.2","MA1234.1","MA1235.1","MA1236.1","MA1238.2","MA1239.1","MA1240.1","MA1241.1","MA1242.1","MA1243.1","MA1244.1","MA1245.2","MA1246.1","MA1247.1","MA1248.1","MA1250.1","MA1251.1","MA1252.1","MA1253.1","MA1256.1","MA1257.1","MA1258.1","MA1259.1","MA1260.1","MA1261.1","MA1262.1","MA1263.1","MA1264.1","MA1265.2","MA1266.1","MA1267.1","MA1268.1","MA1269.1","MA1270.1","MA1272.2","MA1273.1","MA1274.1","MA1275.1","MA1276.1","MA1277.1","MA1278.1","MA1279.1","MA1280.1","MA1281.1","MA1282.1","MA1283.1","MA1284.1","MA1285.1","MA1286.1","MA1287.1","MA1288.1","MA1289.1","MA1290.1","MA1291.1","MA1292.1","MA1293.1","MA1294.1","MA1295.1","MA1296.1","MA1297.1","MA1298.1","MA1299.1","MA1300.1","MA1301.1","MA1302.1","MA1303.1","MA1304.1","MA1305.1","MA1306.1","MA1307.2","MA1308.1","MA1309.1","MA131","MA1310.1","MA1311.2","MA1312.1","MA1313.1","MA1314.1","MA1315.1","MA1316.1","MA1317.1","MA1318.1","MA1320.1","MA1321.1","MA1322.2","MA1323.1","MA1324.1","MA1325.1","MA1326.1","MA1327.2","MA1328.1","MA1329.2","MA1330.1","MA1331.1","MA1332.1","MA1333.1","MA1334.1","MA1335.1","MA1336.1","MA1337.1","MA1338.2","MA1339.1","MA1340.1","MA1341.1","MA1343.1","MA1344.1","MA1345.1","MA1346.1","MA1348.1","MA1349.1","MA1350.1","MA1351.2","MA1352.1","MA1353.1","MA1355.1","MA1356.1","MA1357.1","MA1358.1","MA1359.2","MA1360.2","MA1361.1","MA1362.1","MA1363.1","MA1364.1","MA1365.2","MA1366.1","MA1367.1","MA1368.2","MA1369.1","MA1370.1","MA1371.1","MA1372.1","MA1373.1","MA1374.1","MA1375.1","MA1376.1","MA1377.1","MA1378.1","MA1379.1","MA1380.1","MA1381.1","MA1382.1","MA1383.1","MA1384.1","MA1385.1","MA1386.2","MA1387.1","MA1388.1","MA1389.1","MA1390.2","MA1391.2","MA1392.2","MA1393.2","MA1394.2","MA1396.1","MA1397.1","MA1398.2","MA1400.1","MA1401.1","MA1402.1","MA1403.1","MA1404.1","MA1405.1","MA1406.1","MA1408.1","MA1409.1","MA1410.1","MA1411.1","MA1412.1","MA1414.1","MA1415.1","MA1416.1","MA1417.1","MA1424.1","MA1425.1","MA1426.1","MA1427.2","MA1428.1","MA1430.1","MA1659.1","MA1660.1","MA1661.1","MA1662.1","MA1663.2","MA1664.2","MA1665.2","MA1666.2","MA1667.2","MA1669.1","MA1670.1","MA1671.1","MA1672.1","MA1673.1","MA1674.2","MA1675.1","MA1676.2","MA1677.1","MA1678.2","MA1679.1","MA1681.1","MA1682.1","MA1685.1","MA1686.1","MA1687.1","MA1688.1","MA1689.1","MA1690.1","MA1691.1","MA1692.1","MA1693.1","MA1694.1","MA1695.1","MA1696.1","MA1697.1","MA1698.1","MA1706.1","MA1732.1","MA1733.1","MA1734.1","MA1735.1","MA1736.1","MA1737.1","MA1738.1","MA1739.1","MA1740.1","MA1741.1","MA1742.1","MA1743.1","MA1744.1","MA1745.1","MA1746.1","MA1747.1","MA1748.1","MA1749.1","MA1750.1","MA1751.1","MA1752.1","MA1753.1","MA1754.1","MA1755.1","MA1756.1","MA1757.1","MA1758.1","MA1759.1","MA1760.1","MA1761.1","MA1762.1","MA1764.1","MA1765.1","MA1766.1","MA1767.1","MA1768.1","MA1769.1","MA1770.1","MA1771.1","MA1772.1","MA1773.1","MA1774.1","MA1775.1","MA1776.1","MA1777.1","MA1778.1","MA1779.1","MA1780.1","MA1781.1","MA1783.1","MA1784.1","MA1785.1","MA1786.1","MA1787.1","MA1788.1","MA1789.1","MA1790.1","MA1791.1","MA1792.1","MA1793.1","MA1794.1","MA1795.1","MA1796.1","MA1797.1","MA1799.1","MA1800.1","MA1801.1","MA1802.1","MA1803.1","MA1804.1","MA1805.1","MA1807.1","MA1808.1","MA1809.1","MA1810.1","MA1811.1","MA1812.1","MA1813.1","MA1814.1","MA1815.1","MA1816.1","MA1817.1","MA1818.1","MA1819.1","MA1820.1","MA1821.1","MA1822.1","MA1823.1","MA1824.1","MA1825.1","MA1826.1","MA1827.1","MA1828.1","MA1829.1","MA1830.1","MA1831.1","MA1832.1","MA1833.1","MA1834.1","MA1835.1","MA2005.1","MA2006.1","MA2007.1","MA2008.1","MA2009.1","MA2010.1","MA2011.1","MA2012.1","MA2013.1","MA2014.1","MA2015.1","MA2016.1","MA2017.1","MA2018.1","MA2019.1","MA2020.1","MA2021.1","MA2022.1","MA2023.1","MA2024.1","MA2025.1","MA2026.1","MA2027.1","MA2028.1","MA2029.1","MA2030.1","MA2031.1","MA2032.1","MA2033.1","MA2034.1","MA2035.1","MA2036.1","MA2037.1","MA2038.1","MA2039.1","MA2040.1","MA2041.1","MA2042.1","MA2043.1","MA2044.1","MA2045.1","MA2046.1","MA2047.1","MA2048.1","MA2049.1","MA2050.1","MA2051.1","MA2052.1","MA2053.1"
};

static String[] strand_list = {"plus","minus"};

//pal
static String[] file_name = {"arabidopsis_thaliana.final.annotatePeak.jaspar.AT2G37040.4000up.bed.chip_hub.2.bed",
"arabidopsis_thaliana.final.annotatePeak.jaspar.AT3G10340.4000up.chip_hub.2.bed",
"arabidopsis_thaliana.final.annotatePeak.jaspar.AT3G53260.4000up.chip_hub.2.tab.bed",
"arabidopsis_thaliana.final.annotatePeak.jaspar.AT5G04230.4000up.bed.chip_hub.2.bed"
};

static String[] pal = {"AT2G37040","AT3G10340","AT3G53260","AT5G04230"};

static LinkedHashMap<String, Integer[][][][][]> gene_to_region_motif;
static LinkedHashMap<String, Integer[][][][][]> gene_to_region_tf;

void init_ocr_anal() throws IOException{

	String[] test_file_name = {"AT1G02230.sort.uniq", "AT1G01780.sort.uniq"};
	gene_to_region_motif = new LinkedHashMap();
	gene_to_region_tf = new LinkedHashMap();

	Pattern colon_pattern = Pattern.compile(":");

	Pattern tab_pattern = Pattern.compile("\\t");

	Pattern space_pattern = Pattern.compile("\\s+");

	
	//LinkedList<String> ocr_key = new LinkedList<String>(peak_start.keySet());

	for(int z=0; z<test_file_name.length; z++){

		String id = test_file_name[z];
		id = id.replaceFirst(".sort.uniq","");
		LinkedList<String> temp_start = peak_start.get(id);

		//file is too large. the input file was created with bedtools intersect file 1 file 2.
		//file 1. big/arabidopsis_thaliana.final.annotatePeak. this file can be downloaded from ChIP-HUB
		//file 2. JASPAR2022_araTha1.bb (downloaded from jaspar, converted to bed format using program such as bigbedtobed from ucsc. then because the file was too big, it was splited according to chr number. each chr file was intersected with file 1.
		//in this folder, only a small number of example file are provided
		
		
		String path = "data/chip_hub_jaspar_inter/" + test_file_name[z];
		
		//System.out.println("ocr error:" + test_file_name[z]);

		//value per file (id)
		Integer[][][][][] motif_val = new Integer[type.length][gene_region.length][temp_start.size()][motif_list.length][2];

		Integer[][][][][] tf_val = new Integer[type.length][gene_region.length][temp_start.size()][tf_list.length][2];

		for(int w = 0;w < motif_val.length; w++){

			for(int x = 0;x < motif_val[w].length; x++){

				for(int y = 0;y < motif_val[w][x].length; y++){
		
					for(int t = 0;t < motif_val[w][x][y].length; t++){
					
						for(int e = 0;e < motif_val[w][x][y][t].length; e++){

							motif_val[w][x][y][t][e] = new Integer(0);
						}
					}
				}
			}
		}

		for(int w = 0;w < tf_val.length; w++){

			for(int x = 0;x < tf_val[w].length; x++){

				for(int y = 0;y < tf_val[w][x].length; y++){
		
					for(int t = 0;t < tf_val[w][x][y].length; t++){
					
						for(int e = 0;e < tf_val[w][x][y][t].length; e++){

							tf_val[w][x][y][t][e] = new Integer(0);
						}
					}
				}
			}
		}

		LinkedList<String> start_list = peak_start.get(id);

		List<String> lines  = FileUtils.readLines(new File(path));

		for(int q= 0; q < lines.size(); q++){

			//1	67129	69473	AT1G01140	Intron	Promoter	MA1284.1	+	TCP1	30

			String one_line = lines.get(q).trim();
			String[] tab_split = tab_pattern.split(one_line);
			if(tab_split.length == 10){

				String cur_region = tab_split[4].trim();
				int region_int = -1;
				for(int e = 0; e < gene_region.length; e++){
					if(cur_region.equals(gene_region[e])){

						region_int = e;
						break;

					}
				}

				String cur_type = tab_split[5].trim();
				int type_int = -1;
				for(int e = 0; e < type.length; e++){
					if(cur_type.equals(type[e])){

						type_int = e;
						break;

					}
				}
				
				String cur_start = tab_split[1].trim();
				int start_int = -1;
				for(int e = 0; e < start_list.size(); e++){
					if(cur_start.equals(start_list.get(e))){

						start_int = e;
						break;

					}
				}
				
				String cur_motif = tab_split[6].trim();
				int motif_int = -1;
				for(int e = 0; e < motif_list.length; e++){
					if(cur_motif.equals(motif_list[e])){

						motif_int = e;
						
						break;

					}
				}
				String cur_strand = tab_split[7].trim();

				int strand_int = -1;
				if(cur_strand.equals("-")){

					strand_int = 1;

				}else{

					strand_int = 0;
				}

				String cur_tf = tab_split[8].trim();
				int tf_int = -1;
				for(int e = 0; e < tf_list.length; e++){
					if(cur_tf.equals(tf_list[e])){

						tf_int = e;
						break;

					}
				}
				

				if(type_int != -1 && region_int != -1 && start_int != -1 && motif_int != -1 && strand_int != -1){

					motif_val[type_int][region_int][start_int][motif_int][strand_int] = motif_val[type_int][region_int][start_int][motif_int][strand_int]+1;
				}

				if(type_int != -1 && region_int != -1 && tf_int != -1 && start_int != -1 && strand_int != -1){

					tf_val[type_int][region_int][start_int][tf_int][strand_int] = tf_val[type_int][region_int][start_int][tf_int][strand_int]+1;
				}

				

				
			}//length == 10
		}//q

		gene_to_region_motif.put(id, motif_val);
		gene_to_region_tf.put(id, tf_val);
	}//z
 

}//method


static String[] state ={"S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S26","S27","S28","S29","S30","S31","S32","S33","S34","S35","S36"}; 

static LinkedHashMap<String, String> acc_to_gene_name;

static LinkedHashMap<String, Integer[][]> gene_to_region_state;

void map_gene_pcsd_region_state() throws IOException{

	gene_to_region_state = new LinkedHashMap();
	acc_to_gene_name = new LinkedHashMap();

	Pattern comma_pattern = Pattern.compile(",");

	Pattern tab_pattern = Pattern.compile("\\t");

	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/chrom_gene_state.all.sort";
	
	String[] gene_region_temp = {"downstream","exon","five_prime_UTR","intron","promoter","pseudogenic_exon","three_prime_UTR"};

	List<String> lines  = FileUtils.readLines(new File(path));

	//1	AT1G01260	JAM2	protein_coding	basic helix-loop-helix (bHLH) DNA-binding superfamily protein	chr1	108946	111609	+	exon,three_prime_UTR

	for(int q= 0; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();
		String[] tab_split = tab_pattern.split(one_line);
		String id = tab_split[1].trim();
		String gene_name = tab_split[2].trim();
		acc_to_gene_name.put(id,gene_name);
		
		int cur_state = new Integer(tab_split[0].trim()).intValue()-1;
		Integer[][] freq = new Integer[gene_region_temp.length][state.length];
		
		for(int x = 0; x < freq.length; x++){
			for(int y = 0; y < freq[x].length; y++){
				freq[x][y] = new Integer(0);
			}
		}
				
		String cur_gene_region = tab_split[9].trim();
		String[] comma_split = comma_pattern.split(cur_gene_region);

		for(int x = 0; x < comma_split.length; x++){

			String region_temp = comma_split[x].trim();
			for(int y = 0; y < gene_region_temp.length; y++){

				if(region_temp.equals(gene_region_temp[y])){
					freq[y][cur_state] = new Integer(1);
				}
			}

		}//x

		if(gene_to_region_state.containsKey(id)){

			Integer[][] exist = gene_to_region_state.get(id);

			for(int g  = 0; g < exist.length; g++){

				if(freq[g][cur_state] != 0){
					exist[g][cur_state] = exist[g][cur_state] + freq[g][cur_state];
				}
			}
								
			gene_to_region_state.put(id,exist);
		}else{
			gene_to_region_state.put(id,freq);
								
		}


	}//q

 

}//method

static LinkedHashMap<String,String> sra_indi_to_exp;
void map_sra_exp() throws IOException{

	System.out.println("map_sra_exp():" );
	sra_indi_to_exp = new LinkedHashMap();


 	
  	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	
	String input="data/dhs.exp.table.ara";
	
	List<String> lines  = FileUtils.readLines(new File(input));
	
	for(int q= 0; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();
		
		String[] one_line_space = space_pattern.split(one_line);
			
		String exp = one_line_space[0];
		String indi = one_line_space[1];
			
		sra_indi_to_exp.put(indi,exp);
			
			

	}//q
System.out.println("sra_indi_to_exp:" +sra_indi_to_exp.size());
}

static String[] sra_indi = {
"DRX122176","DRX122177","DRX122178","DRX122179","ERX3444754","ERX3444755","ERX3444756","ERX3444757","ERX3444758","ERX3444759","ERX3444760","ERX3444761","ERX3444762","ERX3444763","ERX3444764","ERX3444765","ERX3444766","ERX3444767","ERX3444768","ERX3444769","ERX3444770","ERX3444771","ERX3444797","ERX3444798","ERX3444799","ERX3444800","ERX3444801","ERX3444802","ERX3444803","ERX3444804","ERX3444805","ERX3444806","ERX3444807","ERX3444808","SRX2528906_SRR5218777","SRX2528907_SRR5218778","SRX2528907_SRR5219624","SRX10089717","SRX10089718","SRX10089719","SRX10089720","SRX10089721","SRX10089722","SRX10089723","SRX10089724","SRX10771698","SRX1096548","SRX1096549","SRX1096550","SRX1096551","SRX1098134","SRX1098135","SRX1098136","SRX1098137","SRX1098138","SRX111004","SRX111005","SRX111006","SRX111007","SRX111008","SRX111009","SRX111010","SRX111011","SRX1204325","SRX1204335","SRX2000799","SRX2000800","SRX2000801","SRX2000802","SRX2000803","SRX2000804","SRX2000805","SRX2000806","SRX2000807","SRX2000808","SRX2000809","SRX2000810","SRX2000812","SRX2311144","SRX2311146","SRX2528905","SRX2528907","SRX2528908","SRX2528909","SRX2528910","SRX277579","SRX277580","SRX277581","SRX277582","SRX277583","SRX3006634","SRX3006635","SRX3006636","SRX3006637","SRX3006644","SRX3006645","SRX3006646","SRX3006647","SRX3006648","SRX3040866","SRX3040867","SRX3040868","SRX3040869","SRX3040870","SRX3041697","SRX3041698","SRX3041699","SRX3041700","SRX3041701","SRX3041702","SRX3348137","SRX3348138","SRX3348139","SRX3348140","SRX3348141","SRX3348142","SRX3348143","SRX3348144","SRX3348145","SRX3348146","SRX3348147","SRX3348148","SRX3348149","SRX3348150","SRX3348151","SRX3348152","SRX3503847","SRX3503848","SRX3503849","SRX3503850","SRX391959","SRX391960","SRX391961","SRX391962","SRX391963","SRX391964","SRX391965","SRX391966","SRX391967","SRX391968","SRX391969","SRX391970","SRX391971","SRX391972","SRX391986","SRX391987","SRX391988","SRX391989","SRX391990","SRX391991","SRX391992","SRX391993","SRX391994","SRX391995","SRX391996","SRX391997","SRX4101257","SRX4101258","SRX4101259","SRX4101260","SRX4101261","SRX4101262","SRX4101263","SRX4382140","SRX4382141","SRX4382142","SRX4382143","SRX4382144","SRX4382145","SRX4382146","SRX4382147","SRX4916677","SRX4916678","SRX4916679","SRX4916680","SRX4916683","SRX4916684","SRX5000426","SRX5000427","SRX5000428","SRX5000429","SRX5000430","SRX5000431","SRX5000432","SRX5000433","SRX5000434","SRX5000435","SRX5000436","SRX5000437","SRX5000438","SRX5000439","SRX5000440","SRX5000441","SRX5000442","SRX5000443","SRX5000444","SRX5000445","SRX5000446","SRX5000447","SRX5000448","SRX5000449","SRX5000450","SRX5000451","SRX5000452","SRX5036313","SRX5036314","SRX5036315","SRX5052456","SRX5052457","SRX5052458","SRX5052459","SRX5052460","SRX5052461","SRX5052462","SRX5052463","SRX5088430","SRX5088431","SRX5088432","SRX5088433","SRX5413160","SRX5413161","SRX5413168","SRX5413169","SRX5534550","SRX5534551","SRX5534552","SRX6081143","SRX6081144","SRX6784996","SRX6784997","SRX6784998","SRX6785001","SRX6785002","SRX6785003","SRX6785004","SRX6785005","SRX6785006","SRX6785007","SRX6785008","SRX6785020","SRX7442790","SRX7442791","SRX7442792","SRX7442793","SRX7442794","SRX7442795","SRX7442796","SRX7442797","SRX7442798","SRX7780553","SRX7780554","SRX7780555","SRX7780559","SRX7780560","SRX7780565","SRX7780566","SRX7780567","SRX7780572","SRX7780573","SRX7780574","SRX7780579","SRX7780580","SRX7780581","SRX8770516","SRX8770517","SRX8840342","SRX8840343","SRX8843250","SRX8843251","SRX8843252","SRX8843253","SRX8843254","SRX8843255","SRX8843256","SRX8844360","SRX8844361","SRX8844362","SRX8844363","SRX8844364","SRX8844365","SRX8844366","SRX8844367","SRX8844368","SRX8844369","SRX8861288","SRX8861289","SRX8861290","SRX8861291","SRX8861292","SRX8861293","SRX8861294","SRX8861295","SRX8861296","SRX8861297","SRX8861298","SRX8861299","SRX8861300","SRX8861301","SRX8861302","SRX8861303","SRX8861304","SRX8861305","SRX8861306","SRX8861307","SRX8861308","SRX8861309","SRX8861310","SRX8861311","SRX8861312","SRX8861313","SRX8861314","SRX8861315","SRX8861316","SRX8861317","SRX8861318","SRX8861319","SRX8861320","SRX8861321","SRX8861322","SRX8861323","SRX8861324","SRX8861325","SRX8861326","SRX8861327","SRX8861328","SRX8861329","SRX8861330","SRX8861331","SRX8861332","SRX8861333","SRX8861334","SRX894596","SRX895212","SRX895220","SRX895226","SRX895237","SRX895244","SRX895266","SRX895269","SRX895270","SRX895271","SRX895316","SRX895317","SRX895318","SRX895319","SRX895320","SRX895321","SRX895322","SRX895323","SRX895324","SRX895325","SRX895326","SRX895327","SRX895328","SRX895329","SRX895330","SRX895331","SRX895332","SRX895333","SRX895334","SRX895336","SRX895355","SRX895356","SRX895359","SRX895361","SRX895364","SRX895366","SRX9074796","SRX9074818","SRX9074819","SRX9074820","SRX9074821","SRX9074822","SRX9132136","SRX9132137","SRX9132138","SRX9191285","SRX9191286","SRX9770773","SRX9770774","SRX9770775","SRX9770776","SRX9770777","SRX9770778","SRX9770779","SRX9770780","SRX9770781","SRX9770782","SRX9770783","SRX9770784","SRX9770785","SRX9770786","SRX9770787","SRX9770788","SRX9770789","SRX9770790","SRX9770791","SRX9770792","SRX9770793","SRX9770794","SRX9770795","SRX9770796","SRX9770797","SRX9770798","SRX9770799","SRX9770800","SRX9770801","SRX9770802","SRX9770803","SRX9770804","SRX9770805","SRX9770806","SRX9770807","SRX9770808","SRX9770809","SRX9770810","SRX9770811","SRX9770812","SRX9770813","SRX9770814","SRX9770815","SRX9770816","SRX9770817","SRX9770818","SRX9770819","SRX9770820","SRX9770821","SRX9770822","SRX9770823","SRX9770824","SRX9770825","SRX9770826","SRX9770827","SRX9770828","SRX9770829","SRX9770830","SRX9770831","SRX9770832","SRX9770833","SRX9770834","SRX9770835","SRX9770836","SRX9770837","SRX9770838","SRX9770839","SRX9770840","SRX9770841","SRX9770842","SRX9770843","SRX9770844","SRX9770845","SRX9770846","SRX9770847","SRX9770848","SRX9770849","SRX9770850","SRX9770851","SRX9770852","SRX9770853","SRX9770854","SRX9770855","SRX9770856","SRX9770857","SRX9770858","SRX9770859","SRX9770860","SRX9770862","SRX9770863","SRX9770864","SRX9770865","SRX9770866","SRX9770867","SRX9770868"
};


static LinkedHashMap<String,Boolean[][][][]> gene_to_sra;
static LinkedHashMap<String,Range<Integer>[][][]> gene_to_range;

void map_gene_sra_exp() throws IOException {
System.out.println("map_gene_sra_exp():" );

	gene_to_sra = new LinkedHashMap();
	gene_to_range = new LinkedHashMap();
	
/*
1	74	729	AT1G01010	Intergenic	Enhancer	SRX8861309;SRX8861300;SRX8861307;SRX2528905;SRX8861299;SRX8861295;SRX8861293;SRX277580;SRX9132137;SRX8861288;SRX8861333;SRX8861334;SRX5088430;SRX8861290;SRX277579;SRX8861303;SRX8843255;SRX391972;SRX8861291;SRX8861319;SRX391987;SRX277581;SRX9132138;SRX8861299;SRX8861293;SRX5088431;SRX8861301;SRX8861297;SRX8861295	
1	74	729	AT1G01010	Intergenic	Enhancer	SRX8861309;SRX8861300;SRX8861307;SRX2528905;SRX8861299;SRX8861295;SRX8861293;SRX277580;SRX9132137;SRX8861288;SRX8861333;SRX8861334;SRX5088430;SRX8861290;SRX277579;SRX8861303;SRX8843255;SRX391972;SRX8861291;SRX8861319;SRX391987;SRX277581;SRX9132138;SRX8861299;SRX8861293;SRX5088431;SRX8861301;SRX8861297;SRX8861295	
*/


	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern semi_colon_pattern = Pattern.compile(";");
	Pattern digit_pattern = Pattern.compile("\\d");

	

	File[] spec_files = new File("data/all_peak").listFiles();
	
	//for(int w = 0; w < pal.length; w++){
	for(int w = 0; w < spec_files.length; w++){
		
		String file_name   = spec_files[w].getName();
		
		String cur_gene =file_name.replaceFirst("arabidopsis_thaliana.final.annotatePeak.","");                                   
		
		//cur_gene =cur_gene.replaceFirst(".bed","");
		
		LinkedList<Range<Integer>> range_list = peak_range.get(cur_gene);
		//System.out.println("error peak:"+cur_gene);

		if(range_list != null){
		
			//System.out.println("error peak:"+cur_gene +":" + range_list.get(0).getMaximum());

			Range<Integer>[][][] range_val = new Range[type.length][gene_region.length][range_list.size()];
			
			
			
			List<String> lines  = FileUtils.readLines(spec_files[w]);
			for(int q= 0; q < lines.size(); q++){
								

				String line = lines.get(q).trim();

				System.out.println(line);
				String[] tab_split = space_pattern.split(line);

				if(tab_split.length == 7){
				
					String cur_region = tab_split[4].trim();
					String cur_type = tab_split[5].trim();
					String cur_start = tab_split[1].trim();
					String cur_end = tab_split[2].trim();
					String sra = tab_split[6].trim();
					
					System.out.println("tab_split.length == 11:"+cur_region +":" + cur_type +":" + cur_start +":" + sra);

					Boolean[][][][] val = new Boolean[type.length][gene_region.length][range_list.size()][sra_indi.length];

					for(int x = 0; x < val.length; x++){

						for(int y = 0; y < val[x].length; y++){

							for(int z = 0; z < val[x][y].length; z++){
							
								for(int e = 0; e < val[x][y][z].length; e++){

									val[x][y][z][e] = new Boolean(false);
								}
							}
						}
					}

					int type_int = -1;
					if(cur_type.equals("Promoter")){

						type_int = 0;
					}else{

						type_int = 1;
					}

					int region_int = -1;

					for(int i = 0; i < gene_region.length; i++){
				
						if(cur_region.equals(gene_region[i])){

							region_int = i;
							break;
						}
					}

					int start_int = -1;
					for(int e = 0; e < range_list.size(); e++){
						Range<Integer> temp_range = range_list.get(e);
						
						Integer start = temp_range.getMinimum();
						Integer end = temp_range.getMaximum();
						
						if(line.contains("\t" + start.toString() + "\t" + end.toString())){

							start_int = e;
							break;

						}
					}
					
					System.out.println("type_int:"+type_int +":" + region_int +":" + start_int);
					
					if(type_int == -1 || region_int == -1 || start_int == -1){
					
					
					}else{
					
						String[] sra_split = semi_colon_pattern.split(sra);

						for(int e = 0; e < sra_split.length; e++){

							String cur_sra = sra_split[e].trim();

							int sra_int = -1;

							for(int i = 0; i < sra_indi.length; i++){
						
								if(cur_sra.equals(sra_indi[i])){

									sra_int = i;

									if(type_int != -1 && region_int != -1){
					
										val[type_int][region_int][start_int][i] = new Boolean(true);
										System.out.println("start_int:"+start_int +":" + sra_int );

									}
									break;
								}
							}

							
						}

						//add range
						range_val[type_int][region_int][start_int] = Range.between(new Integer(cur_start),new Integer(cur_end));
						
						
						gene_to_sra.put(cur_gene,val);
						System.out.println("new:val");
						
						
					}
				}//7
			}//q
			
			gene_to_range.put(cur_gene,range_val);
		}//not  null
	}//w


}

static String[] chip_tissue = {"control","flower","inflorescences","leaf","root","root_hair","root_non-hair","seed","seed_coat","seedling","shoot"};

static String[][] chip_sra_to_tissue;

public void  map_sra_to_tissue() throws IOException{

	Pattern tab_pattern = Pattern.compile("\\t+");

	LinkedList<String>[] temp = new LinkedList[chip_tissue.length];
	chip_sra_to_tissue = new String[chip_tissue.length][];

	for(int i = 0; i < temp.length; i++){

		temp[i] = new LinkedList();

	}

	String path = "data/enhancer.meta";


	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
							
		String line = lines.get(q).trim();
		line.replaceAll("\\s","_");
		String[] tab_split = tab_pattern.split(line);

		String tissue_temp = tab_split[0].trim();
		String sra = tab_split[1].trim();

		
		for(int e = 0; e < chip_tissue.length; e++){
		
			if(tissue_temp.equals(chip_tissue[e])){

				temp[e].add(sra);
				break;		
			}
		}

	}//q
		
	for(int i = 0; i < temp.length; i++){

		if(temp[i] != null){

			chip_sra_to_tissue[i] = new String[temp[i].size()];
			chip_sra_to_tissue[i] = temp[i].toArray(new String[0]);
		}
		
	}


}

static LinkedHashMap<String,Boolean[][][][]> gene_to_tissue;

void map_gene_to_tissue() throws IOException {
	System.out.println("map_map_gene_to_tissue():" );

	gene_to_tissue = new LinkedHashMap();
	
	String header = "";
		
	for(int j = 0; j < type.length; j++){//type

		for(int l = 0; l < gene_region.length; l++){//region
				
			for(int e = 0; e < chip_tissue.length; e++){

				header = header + type[j] + "_" + gene_region[l] + "_" + chip_tissue[e] +  ",";
					

			}
		}
	}
		
	System.out.println("gene_to_tissue_header:" + header);
		
		
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern semi_colon_pattern = Pattern.compile(";");
	Pattern digit_pattern = Pattern.compile("\\d");
	
	List<String> keys = new LinkedList<String>(gene_to_sra.keySet());
	System.out.println("gene_to_sra.keySet():" + keys.size() + ":" + keys.get(0));
	
	for(int i = 0; i < keys.size(); i++){

		String acc = keys.get(i);
		LinkedList<String> start_list = peak_start.get(acc);
		Boolean[][][][] val = gene_to_sra.get(acc);
		Boolean[][][][] new_val = new Boolean[val.length][][][];


		for(int j = 0; j < val.length; j++){//type

			new_val[j] = new Boolean[val[j].length][][];

			for(int l = 0; l < val[j].length; l++){//region
				
				new_val[j][l] = new Boolean[start_list.size()][];

				for(int k = 0; k < val[j][l].length; k++){//start
				
					new_val[j][l][k] = new Boolean[chip_tissue.length];
				
					for(int e = 0; e < chip_tissue.length; e++){

						new_val[j][l][k][e] = new Boolean(false);
					}

				}
			}
		}
		

		for(int j = 0; j < val.length; j++){	

			for(int k = 0; k < val[j].length; k++){

				for(int l = 0; l < val[j][k].length; l++){
				
					for(int e = 0; e < val[j][k][l].length; e++){

						if(val[j][k][l][e].booleanValue()){

							String cur_sra_exp = sra_indi[e];

							for(int m = 0; m < chip_sra_to_tissue.length; m++){

								if(chip_sra_to_tissue[m] != null){

									for(int n = 0; n < chip_sra_to_tissue[m].length; n++){
			
											if(cur_sra_exp.equals(chip_sra_to_tissue[m][n])){

											new_val[j][k][l][m] = new Boolean(true);

											break;
										}
									}//n
								}
								
								if(new_val[j][k][l][m]){

									break;
								}
							}//m
						}//if
					}//e
				}//l
			}//k
		}//j	

		gene_to_tissue.put(acc,new_val);
	}//i key
	
	keys = new LinkedList<String>(gene_to_tissue.keySet());
	for(int i = 0; i < keys.size(); i++){

		String acc = keys.get(i);
		LinkedList<String> start_list = peak_start.get(acc);
		Boolean[][][][] temp_tissue = gene_to_tissue.get(acc);
		for(int e = 0; e < type.length; e++){

			for(int q = 0; q < gene_region.length; q++){
			
				for(int w = 0; w < start_list.size(); w++){
					
					for(int r = 0; r < chip_tissue.length; r++){
					
						System.out.println("gene_to_tissue:" + type[e] + ":" + gene_region[q] + ":" + start_list.get(w) + ":" + chip_tissue[r] + ":" + temp_tissue[e][q][w][r].booleanValue());

					}
				}
			}
		}
	}


}//method

///

String[] study_16_pheno ={"As75","Cd111","Co59","Cu65","K39","Mg25","Mn55","Mo98","Ni60","P31","Rb85","S34","Se82","Zn66"};

String[] epi = {"CG","CHG","CHH"};


//get ewas_matal_pheno_per_promoter_enhancer

LinkedHashMap<String,Integer[][][][][]> gene_to_metal_pheno_freq;
void map_ewas_metal_pheno_promoter_enhancer() throws IOException {

	System.out.println("map_ewas_metal_pheno_promoter_enhancer():" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern semi_colon_pattern = Pattern.compile(";");
	Pattern digit_pattern = Pattern.compile("\\d");
	
	gene_to_metal_pheno_freq = new LinkedHashMap();
	
	File[] spec_files = new File("big/ewas/indi2").listFiles();
	
	for(int w = 0; w < spec_files.length; w++){
		
		String file_name   = spec_files[w].getName();
		//ewas.gff.updown3000.bed.AT1G47890.bed
		String cur_gene =file_name.replaceFirst("ewas.gff.updown3000.bed.","");                                   
		
		cur_gene =cur_gene.replaceFirst(".bed","");
		
		Range<Integer>[][][] gene_range = gene_to_range.get(cur_gene);
		//System.out.println("error peak:"+cur_gene);

		if(gene_range != null){
		
			//System.out.println("error peak:"+cur_gene +":" );

			Integer[][][][][] freq_val = new Integer[type.length][gene_region.length][gene_range[0][0].length][study_16_pheno.length][epi.length];
			
			for(int x = 0; x < freq_val.length; x++){

				for(int y = 0; y < freq_val[x].length; y++){

					for(int z = 0; z < freq_val[x][y].length; z++){
							
						for(int e = 0; e < freq_val[x][y][z].length; e++){
							for(int f = 0; f < freq_val[x][y][z][e].length; f++){
								freq_val[x][y][z][e][f] = new Integer(0);
								
							}
						}
					}
				}
			}
					
			
			List<String> lines  = FileUtils.readLines(spec_files[w]);
			for(int q= 0; q < lines.size(); q++){
								
				/*
				1	640	641	.	.	.	CG	Mo98	-0.458691009318133	7.64176900530426e-34	1.17674064976704e-32	study16	-0.452591986946152	6.88377227097905e-33	9.30338350779742e-32	molybdenum concentration	1	630	8899	.	+	.	gene	UPDOWN_ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010	1
1	640	641	.	.	.	CG	Na23	0.519022482315654	4.66765481962001e-44	1.8511367064696e-43	study16	0.426908491983456	7.45515512769258e-29	1.62047218050285e-28	sodium concentration	1	630	8899	.	+	.	gene	UPDOWN_ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010	1
*/

				String line = lines.get(q).trim();

				System.out.println(line);
				String[] tab_split = space_pattern.split(line);

				if(tab_split.length == 25){
				
					String cur_start = tab_split[1].trim();
					String cur_end = tab_split[2].trim();
					String cur_epi = tab_split[6].trim();
					String cur_pheno = tab_split[7].trim();
					
					System.out.println("tab_split.length == 25:"+cur_start +":" + cur_end +":" + cur_epi +":" + cur_pheno);

					
					int epi_int = -1;

					for(int i = 0; i < epi.length; i++){
				
						if(cur_epi.equals(epi[i])){

							epi_int = i;
							break;
						}
					}
					
					int pheno_int = -1;

					for(int i = 0; i < study_16_pheno.length; i++){
				
						if(cur_pheno.equals(study_16_pheno[i])){

							pheno_int = i;
							break;
						}
					}

					System.out.println("epi_int:"+epi_int +":" + pheno_int);
					
					if(epi_int == -1 || pheno_int == -1){
					
					
					}else{
					
						Range<Integer> cur_range = Range.between(new Integer(cur_start), new Integer(cur_end));
					
						for(int f = 0; f < gene_range.length; f++){
						
							boolean find = false;
							for(int g = 0; g < gene_range[f].length; g++){
							
								for(int d = 0; d < gene_range[f][g].length; d++){
								
									if(gene_range[f][g][d] != null){
										if(gene_range[f][g][d].isOverlappedBy(cur_range)){
				
											freq_val[f][g][d][pheno_int][epi_int] = freq_val[f][g][d][pheno_int][epi_int]+1;
											find = true;
											break;
										
									
										}
									}
								}//d
								
								if(find){
								
									break;
								}
							}//region g
						}//type f
				
						
					}
				}//7
			}//q
			
			gene_to_metal_pheno_freq.put(cur_gene,freq_val);
		}//not  null
	}//w
	
}//method




static LinkedHashMap<String,LinkedList<String>> motif_acc_to_domain;
static LinkedHashMap<String,LinkedList<String>> domain_to_motif_acc;

static LinkedHashMap<String,Boolean[]> motif_acc_to_domain_bool;
static LinkedHashMap<String,Boolean[]> domain_to_motif_acc_bool;

static LinkedHashMap<String,Boolean[]> motif_acc_to_cluster;
static LinkedHashMap<String,Boolean[]> domain_to_cluster;

static LinkedList<String>[] cluster_memb_id;
static LinkedList<String>[] cluster_memb_acc;


void  map_motif_acc_to_domain() throws IOException{

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern tri_colon_pattern = Pattern.compile(":::");
	Pattern left_pattern = Pattern.compile("left");
	Pattern semi_pattern = Pattern.compile(";");

	motif_acc_to_domain = new LinkedHashMap();
	domain_to_motif_acc = new LinkedHashMap();
	motif_acc_to_cluster = new LinkedHashMap();
	domain_to_cluster = new LinkedHashMap();
	motif_acc_to_domain_bool = new LinkedHashMap();
	domain_to_motif_acc_bool = new LinkedHashMap();
	cluster_memb_id = new LinkedList[47];
	cluster_memb_acc = new LinkedList[47];

	for(int k = 0; k < 47; k++){

		cluster_memb_id[k] = new LinkedList();
		cluster_memb_acc[k] = new LinkedList();
	}

	String path = "data/JASPAR_2022_matrix_clustering_plants_CORE_cluster_root_motifs.tf";

	List<String> lines  = FileUtils.readLines(new File(path));

	//AC  cluster_10
					
	for(int k = 0; k < lines.size(); k++){
			
		String one_line = lines.get(k).trim();

		if(one_line.contains("AC  cluster_")){

			//System.out.println("one_line:" + one_line);

			int clu_num = new Integer(one_line.replaceFirst("AC  cluster_","").trim()).intValue() -1;
			//System.out.println("clu_num:" + clu_num);

			LinkedList<String> acs = new LinkedList();

			for(int j = k+1; j < lines.size(); j++){
			
				String next_line = lines.get(j).trim();

				if(next_line.contains("merged_AC:")){
				//merged_AC: 'JASPAR_2022_plants_CORE:::MA1284_1,JASPAR_2022_plants_CORE:::MA1062_2,JASPAR_2022_plants_CORE:::MA1065_2,JASPAR_2022_plants_CORE:::MA1283_1,JASPAR_2022_plants_CORE:::MA1285_1'

					//System.out.println("merged_AC:");

					if(next_line.contains(",")){

						//System.out.println("merged_AC:,");

						String[] comma_split = comma_pattern.split(next_line);

						for(int i = 0; i < comma_split.length; i++){

							//System.out.println("81:" + comma_split[i].trim());

							String[] tri_colon_split = tri_colon_pattern.split(comma_split[i].trim());

							String input_motif = tri_colon_split[1].trim();
							input_motif = input_motif.replaceFirst("_",".");

							acs.add(input_motif);

						}

					}else{

						String[] tri_colon_split = tri_colon_pattern.split(next_line);
						acs.add(tri_colon_split[1].trim());
						//System.out.println("93:" + tri_colon_split[0].trim());
						
					}
				

					break;
				}//ac
			}//j

			LinkedList<String> ids = new LinkedList();

			for(int j = k+1; j < lines.size(); j++){
			
				String next_line = lines.get(j).trim();

				if(next_line.contains("merged_ID:")){

					//System.out.println("merged_ID:");
				//merged_AC: 'JASPAR_2022_plants_CORE:::MA1284_1,JASPAR_2022_plants_CORE:::MA1062_2,JASPAR_2022_plants_CORE:::MA1065_2,JASPAR_2022_plants_CORE:::MA1283_1,JASPAR_2022_plants_CORE:::MA1285_1'
					next_line = next_line.replaceFirst("CC  merged_ID: ","").trim();
					next_line = next_line.replaceAll("'","");
					if(next_line.contains(",")){

						//System.out.println(",");

						String[] comma_split = comma_pattern.split(next_line);

						for(int i = 0; i < comma_split.length; i++){

							//System.out.println("122:" +comma_split[0].trim() );

							ids.add(comma_split[i].trim());

						}

					}else{

						ids.add(next_line);
						//System.out.println("131:" +next_line);
						
					}
					break;
				}//id
			}//j

			cluster_memb_id[clu_num].addAll(ids);
			cluster_memb_acc[clu_num].addAll(acs);

			for(int j = 0; j < acs.size();j++){

				String acc = acs.get(j);
				String id = ids.get(j);
				//System.out.println("145:" +acc + ":" + id);

				if(motif_acc_to_domain.containsKey(acc)){

					LinkedList<String> exist = motif_acc_to_domain.get(acc);
					if(!exist.contains(id)){
						//System.out.println("error duplicate acc:" + acc + ":" + id);
						exist.add(id);
					}
					motif_acc_to_domain.put(acc,exist);
				}else{

					LinkedList<String> exist = new LinkedList();
					exist.add(id);
					motif_acc_to_domain.put(acc,exist);
					//System.out.println("145:" +acc + ":" + id);
				}

				if(domain_to_motif_acc.containsKey(id)){

					LinkedList<String> exist = domain_to_motif_acc.get(id);
					if(!exist.contains(acc)){
						exist.add(acc);
					}
					domain_to_motif_acc.put(id,exist);
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(acc);
					domain_to_motif_acc.put(id,exist);
				}
			}//j
		}
	}

	

	for(int k = 0; k < 47; k++){

		for(int i = 0; i < cluster_memb_id[k].size(); i++){


			String cur_id = cluster_memb_id[k].get(i);
			int clus_num = k;
			if(domain_to_cluster.containsKey(cur_id)){

				Boolean[] exist = domain_to_cluster.get(cur_id);
				exist[clus_num] = new Boolean(true);
				domain_to_cluster.put(cur_id,exist);
			}else{

					
				Boolean[] exist = new Boolean[47];
				for(int g = 0; g < exist.length; g++){

					exist[g] = new Boolean(false);
				}
				exist[clus_num] = new Boolean(true);
				domain_to_cluster.put(cur_id,exist);
			}
		}
	}
	

	
	for(int k = 0; k < 47; k++){
		

		System.out.println("cluster_memb_acc:cluster"+(k+1)+":"+Arrays.toString(cluster_memb_acc[k].toArray()));

		for(int i = 0; i < cluster_memb_acc[k].size(); i++){


			String cur_acc = cluster_memb_acc[k].get(i);
			int clus_num = k;
			if(motif_acc_to_cluster.containsKey(cur_acc)){

				Boolean[] exist = motif_acc_to_cluster.get(cur_acc);
				exist[clus_num] = new Boolean(true);
				motif_acc_to_cluster.put(cur_acc,exist);
			}else{

					
				Boolean[] exist = new Boolean[47];
				for(int g = 0; g < exist.length; g++){

					exist[g] = new Boolean(false);
				}
				exist[clus_num] = new Boolean(true);
				motif_acc_to_cluster.put(cur_acc,exist);
			}
			
		}
	}

	LinkedList<String> all_id = new LinkedList();
				
	
	for(int k = 0; k < 47; k++){

		System.out.println("cluster_memb_id:cluster"+(k+1)+":"+Arrays.toString(cluster_memb_id[k].toArray()));
		all_id.addAll(cluster_memb_id[k]);

	}

	Set<String> pfams_h_inv = new LinkedHashSet<>();
	pfams_h_inv.addAll((List)all_id);
					
	// Clear the list
	((List)all_id).clear();
											  
	// add the elements of set
	// with no duplicates to the list
	((List)all_id).addAll(pfams_h_inv);


	List<String> keys = new LinkedList<String>(motif_acc_to_domain.keySet());
	for(int i = 0; i < keys.size(); i++){

		String acc = keys.get(i);
		LinkedList<String> id = motif_acc_to_domain.get(keys.get(i));
		System.out.println("acc:" + acc + ":" + Arrays.toString(id.toArray()));
		Boolean[] id_pres = new Boolean[all_id.size()];
		for(int q = 0; q < all_id.size(); q++){

			id_pres[q] = new Boolean(false);
		}

		for(int w = 0; w < id.size(); w++){

			String cur_id = id.get(w);

			for(int q = 0; q < all_id.size(); q++){

				String comp_id = all_id.get(q);
				if(cur_id.equals(comp_id)){

					id_pres[q] = new Boolean(true);
					break;
				}
			}
		}//w

		motif_acc_to_domain_bool.put(acc,id_pres);
	}

	LinkedList<String> all_acc = new LinkedList();
				
	
	for(int k = 0; k < 47; k++){

		System.out.println("cluster_memb_acc:cluster"+(k+1)+":"+Arrays.toString(cluster_memb_acc[k].toArray()));
		all_acc.addAll(cluster_memb_acc[k]);

	}

	pfams_h_inv = new LinkedHashSet<>();
	pfams_h_inv.addAll((List)all_acc);
					
	// Clear the list
	((List)all_acc).clear();
											  
	// add the elements of set
	// with no duplicates to the list
	((List)all_acc).addAll(pfams_h_inv);


	keys = new LinkedList<String>(domain_to_motif_acc.keySet());
	for(int i = 0; i < keys.size(); i++){

		String id = keys.get(i);
		LinkedList<String> acs = domain_to_motif_acc.get(keys.get(i));
		System.out.println("id:" + id + ":" + Arrays.toString(acs.toArray()));

		Boolean[] acc_pres = new Boolean[all_acc.size()];

		for(int q = 0; q < all_acc.size(); q++){

			acc_pres[q] = new Boolean(false);
		}



		for(int w = 0; w < acs.size(); w++){

			String cur_acc = acs.get(w);

			for(int q = 0; q < all_acc.size(); q++){

				String comp_acc = all_acc.get(q);
				if(cur_acc.equals(comp_acc)){

					acc_pres[q] = new Boolean(true);
					break;
				}
			}
		}//w

		domain_to_motif_acc_bool.put(id,acc_pres);
	}
	

}

//create atac-seq width (width is end position - start position) in phloem tissue

static String[] ranges = {"100","200","300","400","500","600","700","800","900","1000","1500","2000","g2000"};

static LinkedHashMap<String,Boolean[]> atac_to_width;

void  calc_atac_width_phloem() throws IOException{

	atac_to_width = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	atac_to_width = new LinkedHashMap();
	
	String path = "data/out_bedintersect_a_b";
	
	List<String> lines  = FileUtils.readLines(new File(path));
	//5	26886829	26887655	.	.	.	0.00014053592288904	0.433798923983554	0.147	0.067	1	8	AT5G67385
//5	26889027	26891586	.	.	.	7.11E-07	0.352468319163545	0.328	0.225	0.0161662053073882	7	AT5G67390

	List<String> gff_key = new LinkedList<String>(atac_to_width.keySet());	
	System.out.println("gff_key:" +gff_key.size());
		
	for(int k = 3; k <lines.size(); k++){
				
		String one_line = lines.get(k);

		//System.out.println("one_line:" + one_line);
		
		String[] space_split = tab_pattern.split(one_line);
		
		String chr = space_split[0].trim();
		
		if(space_split[1].trim().contains(".")){
		
			space_split[1] = space_split[1].trim().substring(0,space_split[1].trim().indexOf("."));
		}
		Integer start = new Integer(space_split[1].trim());
		
		if(space_split[2].trim().contains(".")){
		
			space_split[2] = space_split[2].trim().substring(0,space_split[2].trim().indexOf("."));
		}
		
		Integer end = new Integer(space_split[2].trim());
			
		int width = end.intValue() - start.intValue();
		
		Boolean[] value = new Boolean[ranges.length];
		for(int q = 0; q < value.length; q++){
			 
			 value[q] = new Boolean(false);
		}

		for(int q = 0; q < ranges.length; q++){
		
			if(width >= 1 && width < 100){
			
				value[0] = new Boolean(true);
			
			}else if(width >= 100 && width < 200){
			
				value[1] = new Boolean(true);
			
			}else if(width >= 200 && width < 300){
			
				value[2] = new Boolean(true);
			
			}else if(width >= 300 && width < 400){
			
				value[3] = new Boolean(true);
			
			}else if(width >= 400 && width < 500){
			
				value[4] = new Boolean(true);
			
			}else if(width >= 500 && width < 600){
			
				value[5] = new Boolean(true);
			
			}else if(width >= 600 && width < 700){
			
				value[6] = new Boolean(true);
			
			}else if(width >= 700 && width < 800){
			
				value[7] = new Boolean(true);
			
			}else if(width >= 800 && width < 900){
			
				value[8] = new Boolean(true);
			
			}else if(width >= 900 && width < 1000){
			
				value[9] = new Boolean(true);
			
			}else if(width >= 1000 && width < 1500){
			
				value[10] = new Boolean(true);
			
			}else if(width >= 1500 && width < 2000){
			
				value[11] = new Boolean(true);
			
			}else if(width >= 2000){
			
				value[12] = new Boolean(true);
			
			}
			 
			 
		}
		String gene = space_split[space_split.length-1].trim();
		
		if(atac_to_width.containsKey(gene)){
			
			Boolean[] exist = atac_to_width.get(gene);
			for(int r = 0; r < value.length; r++){
			
				if(value[r]){
					exist[r] = true;
				}
			}
			atac_to_width.put(gene,exist);
				
		}else{
			
			atac_to_width.put(gene,value);
		}
		
	}//k
	
	
	
	List<String> keys_list = new LinkedList<String>(atac_to_width.keySet());

	for(int q = 0; q < keys_list.size(); q++){
	
		String key = keys_list.get(q);
		System.out.println("final:"+key + ":" + Arrays.toString(atac_to_width.get(key)));
	}
}//method

static LinkedHashMap<String, LinkedList<String>> gene_to_phloem_atac;

void create_phloem_atac_seq() throws IOException{

	System.out.println("create_phloem_atac_seq:");
	Pattern space_pattern = Pattern.compile("\\s+");
	
	gene_to_phloem_atac = new LinkedHashMap();
	
	String input="data/arabidopsis_atac_phloem";
	
	List<String> lines  = FileUtils.readLines(new File(input));
	/*
	1	0	15760	.	.	.	5.70E-10	-0.267160354408712	0.088	0.158	1.30E-05	0	AT1G01010
1	13483	15760.1	.	.	.	2.64E-10	-0.27266996864857	0.063	0.158	6.02E-06	1	AT1G01030
1	13483	15760.2	.	.	.	2.99E-10	-0.293172886540679	0.071	0.156	6.80E-06	2	AT1G01030
1	13483	15760.3	.	.	.	1.71E-13	0.440797510735105	0.249	0.134	3.90E-09	4	AT1G01030

*/
	
	for(int q= 3; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();
		
		if(!one_line.equals("")){
		
			String[] one_line_tab = space_pattern.split(one_line);
			
			
			String chr = one_line_tab[0].trim();
			
			String start = one_line_tab[1].trim();
			if(start.contains(".") && !start.contains("E")){
			
				start = start.substring(0, start.indexOf("."));
			}
			String end = one_line_tab[2].trim();
			if(end.contains(".") && !start.contains("E")){
			
				end = end.substring(0, end.indexOf("."));
			}
			
			one_line_tab[1] = new Integer(new BigDecimal(start).intValue()).toString();
			one_line_tab[2] = new Integer(new BigDecimal(end).intValue()).toString();
			
			String key = one_line_tab[one_line_tab.length-1].trim();
			
			String val = chr + "\t" + one_line_tab[1] + "\t" + one_line_tab[2];
			
			if(gene_to_phloem_atac.containsKey(key)){

				LinkedList<String> exist = gene_to_phloem_atac.get(key);
				if(!exist.contains(val)){
					exist.add(val);
				}
				gene_to_phloem_atac.put(key,exist);
			}else{
				LinkedList<String> exist = new LinkedList();
				exist.add(val);
				gene_to_phloem_atac.put(key,exist);
			}
			
		}
	}

}//lines




static LinkedHashMap<String,String> name_to_seq;
static LinkedHashMap<String,String> name_to_id;
static LinkedHashMap<String,String> seq_to_id;
static LinkedHashMap<String,String> id_to_name;

void get_cis_reg_list() throws IOException{
        
        Pattern space_pattern = Pattern.compile("\\s+");
	name_to_seq = new LinkedHashMap();
	name_to_id = new LinkedHashMap();
        seq_to_id = new LinkedHashMap();
	id_to_name  = new LinkedHashMap();
         
        String in_path = "data/place_dat.txt";
        List<String> lines  = FileUtils.readLines(new File(in_path));
        
        //UPRMOTIFIIAT         CCNNNNNNNNNNNNCCACG            19 S000426
	//VOZATVPP             GCGTNNNNNNNACGC                15 S000456

	for(int w = 0; w < lines.size(); w++){
            
            String temp = lines.get(w).trim();
            
            if(!temp.equals("")){
                
                //if(temp.startsWith(">")){
                    
                    String[] split_temp = space_pattern.split(temp);
		    String name = split_temp[0].trim();
		    String seq = split_temp[1].trim();
		    String id = split_temp[3].trim();
		    name_to_seq.put(name,seq);
		    name_to_id.put(name,id);
		    seq_to_id.put(seq,id);
		    seq_to_id.put(id,name);
		}
	}//w
                    
}

static LinkedHashMap<String,String> id_to_desc;

void map_id_to_desc() throws IOException{
        
        Pattern space_pattern = Pattern.compile("\\s+");
	id_to_desc = new LinkedHashMap();
	
         
        String in_path = "data/place_seq.txt";
        List<String> lines  = FileUtils.readLines(new File(in_path));
        
        /*

	ID   -10PEHVPSBD
XX
AC   S000392
XX
DT   20-Feb-2002 (last modified) uchi
XX
DE   "-10 promoter element" found in the barley (H.v.) chloroplast
DE   psbD gene promoter; Involved in the expression of the plastid
DE   gene psbD which encodes a photosystem II reaction center
DE   chlorophyll-binding protein that is activated by blue, white or
DE   UV-A light;
XX
KW   psbD; chloroplast gene expression; circadian rhythms; light
KW   regulation;
XX
OS   barley (Hordeum vulgare)
XX
RA   Thum KE, Kim M, Morishige DT, Eibl C, Koop HU, Mullet JE
RT   Analysis of barley chloroplast psbD light-responsive promoter
RT   elements in transplastomic tobacco
RL   Plant Mol Biol  47: 353-366 (2001)
RD   PubMed: 11587507;
XX
SQ
     TATTCT
//
*/
	for(int w = 0; w < lines.size(); w++){
            
            String temp = lines.get(w).trim();
            
            if(temp.startsWith("AC")){

		String id = temp.replaceFirst("AC","").trim();

		String desc = "";
		boolean find = false;
                
                for(int j = w+1; j < lines.size(); j++){
            
            		String next_str = lines.get(j).trim();

			if(next_str.equals("//")){
				break;
			}
                    
		        if(next_str.startsWith("DE")){

				desc = desc + next_str.replaceFirst("DE","").trim() + " ";
				find = true;
			}else{

				if(find){

					break;
				}
			}
		}

		if(desc.equals("")){

			desc = "NA";

		}

		id_to_desc.put(id,desc);
		break;
	    }//ac
	}//w
                    
}

    


static String[] jaspar_array = {"$MA0001.2","$MA0005.2","$MA0008.1","$MA0110.2","$MA0121.1","$MA0548.1","$MA0549.1","$MA0550.1","$MA0551.1","$MA0552.1","$MA0553.1","$MA0554.1","$MA0555.1","$MA0556.1","$MA0557.1","$MA0558.1","$MA0559.1","$MA0560.1","$MA0561.1","$MA0562.1","$MA0563.1","$MA0564.1","$MA0565.1","$MA0566.1","$MA0567.1","$MA0568.1","$MA0569.1","$MA0570.1","$MA0571.1","$MA0572.1","$MA0573.1","$MA0574.1","$MA0575.1","$MA0576.1","$MA0577.1","$MA0578.1","$MA0579.1","$MA0580.1","$MA0581.1","$MA0582.1","$MA0583.1","$MA0584.1","$MA0585.1","$MA0586.1","$MA0587.1","$MA0588.1","$MA0589.1","$MA0590.1"};

static LinkedHashMap<String,Double[][]> place_to_jaspar;

void map_place_seq_to_jaspar() throws IOException{

	Pattern space_pattern = Pattern.compile("\\s+");
	place_to_jaspar = new LinkedHashMap();
	
         
        String in_path = "data/all4.Rout";
        List<String> lines  = FileUtils.readLines(new File(in_path));
        
        /*

[1] "TATTCT"
$MA0008.1
   score relScore 
 8.51200 70.93333 

$MA0121.1
    score  relScore 
 8.213333 68.444443 

*/
	for(int w = 0; w < lines.size(); w++){
            
            String temp = lines.get(w).trim();
            
            if(temp.startsWith("\\[1\\]")){

		String id = temp.replaceFirst("\\[1\\]","").trim();
		System.out.println("699 id:" + id);
		id = id.replaceAll("\"","").trim();
		System.out.println("701 id:" + id);
		Double[][] val = new Double[jaspar_array.length][2];
		for(int r = w+1; r < val.length; r++){

			for(int e = 0; e < 2; e++){

				val[r][e] = new Double(0.0);
			}
		}

		for(int r = w+1; r < lines.size(); r++){

			if(temp.startsWith("\\[1\\]")){

				place_to_jaspar.put(id,val);
				break;
			}
			
			String next = lines.get(r).trim();
			if(next.startsWith("$")){

				String motif = next;
				int motif_ind = -1;

				for(int e = 0; e < jaspar_array.length; e++){

					if(next.equals(jaspar_array[e].trim())){

						motif_ind = e;
						break;
					}
				}

				String score = lines.get(r+2).trim();
				String[] score_split = space_pattern.split(score);
				val[motif_ind][0] = new Double(score_split[0].trim());
				val[motif_ind][1] = new Double(score_split[1].trim());
			}
		}//r
	}//if
	}//w

    }

static String[] chr_seq_array;
    
void get_chr_seq_array() throws IOException{

	String[] header = {">Chr1 CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02",
	   ">Chr2 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02", 
	   ">Chr3 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02", 
	   ">Chr4 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02", 
	   ">Chr5 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02"
	};
	chr_seq_array = new String[header.length]; 
        
	int[] index = {0,385163,634510,931471,1166726,1508190};

	   
       String in_path = "big/TAIR10_chr_all.fas";
      

	List<String> lines  = FileUtils.readLines(new File(in_path));

	for(int i = 0; i < index.length-1; i++){

		for(int w = index[i]+1; w < index[i+1]; w++){

			chr_seq_array[i] = chr_seq_array[i] + lines.get(w).trim();
		}

		chr_seq_array[i] = chr_seq_array[i].replaceAll(" ","");
	}

                    
}//method



static LinkedHashMap<String,LinkedList<Integer>> gene_to_tss_index;
static LinkedHashMap<String,String> gene_to_strand;
static LinkedHashMap<String,LinkedList<String>> gene_to_utr5_seq;
static LinkedHashMap<String,LinkedList<String>> gene_to_gene_body_seq;

void get_promoter_gene_body() throws IOException{
        
        Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern semicolon_pattern = Pattern.compile(";");

	
	gene_to_strand = new LinkedHashMap();
	gene_to_tss_index = new LinkedHashMap();
	gene_to_utr5_seq = new LinkedHashMap();
	gene_to_gene_body_seq = new LinkedHashMap();
         
        String in_path = "big/TAIR10_GFF3_genes.txt";
        List<String> lines  = FileUtils.readLines(new File(in_path));
        

	for(int w = 0; w < lines.size(); w++){
            
            String temp = lines.get(w).trim();
            
            if(!temp.equals("")){
                
                //if(temp.startsWith(">")){
                    
                    String[] split_temp = tab_pattern.split(temp);
		    String type = split_temp[1].trim();

		    String[] id_split = semicolon_pattern.split(split_temp[8].trim());
		    String gene_name = id_split[0].trim();
		    gene_name = gene_name.replaceFirst("ID=","");
		

		    if(type.equals("gene")){

			//get gene index
			String chr = split_temp[0].trim();
		    	chr = chr.replaceFirst("Chr","");
			int chr_int = new Integer(chr).intValue()-1;
			String strand = split_temp[6].trim();
			//map gene to strand
			gene_to_strand.put(gene_name,strand);

			Integer start = -1;
			Integer end = -1;
			int start_int = -1;
			int end_int = -1;
			String gene_seq = "";
			start_int = new Integer(split_temp[3]).intValue() -1;
			start = new Integer(start_int);
			end = new Integer(split_temp[4]);
			end_int = end.intValue();
			gene_seq = chr_seq_array[chr_int].substring(start_int,end_int);

			if(strand.equals("+")){
		
				
				
			}else{
				

				String gene_seq_rc = reverse_complement(gene_seq);
				gene_seq = gene_seq_rc;

				int rc_start_int = gene_seq.length() - end_int;
				int rc_end_int = gene_seq.length() - start_int +1;

				start_int = rc_start_int;
				end_int = rc_end_int;
				start = new Integer(start_int);
				end = new Integer(end_int);
			}

			for(int j = w+1; j < lines.size(); j++){
            
			    String next = lines.get(w).trim();
			    String[] split_next = tab_pattern.split(next);
		    	    String next_type = split_next[1].trim();

			    

			    if(temp.equals("gene")){

				break;
			    }else if(next_type.equals("five_prime_UTR")){

				Integer tss = new Integer(split_next[4]);
				int tss_int = tss.intValue();
				String utr5_seq = gene_seq.substring(start,tss_int);
				String gene_body_seq = gene_seq.substring(tss_int,gene_seq.length());

				if(strand.equals("+")){

					
						
					
				
				}else{
					tss = gene_seq.length() -(new Integer(split_next[3])-2) -1;
					tss_int = tss.intValue();
					//map gene to tss
					

					utr5_seq = gene_seq.substring(start,tss_int);
					gene_body_seq = gene_seq.substring(tss_int,gene_seq.length());
				}

				if(gene_to_tss_index.containsKey(gene_name)){

					LinkedList<Integer> exist = gene_to_tss_index.get(gene_name);
					if(!exist.contains(tss)){
						exist.add(tss);
					}
					gene_to_tss_index.put(gene_name,exist);
				}else{
					LinkedList<Integer> exist = new LinkedList();
					exist.add(tss);
					gene_to_tss_index.put(gene_name,exist);
				}

				if(gene_to_utr5_seq.containsKey(gene_name)){

					LinkedList<String> exist = gene_to_utr5_seq.get(gene_name);
					if(!exist.contains(utr5_seq)){
						exist.add(utr5_seq);
					}
					gene_to_utr5_seq.put(gene_name,exist);
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(utr5_seq);
					gene_to_utr5_seq.put(gene_name,exist);
				}

				if(gene_to_gene_body_seq.containsKey(gene_name)){

					LinkedList<String> exist = gene_to_gene_body_seq.get(gene_name);
					if(!exist.contains(gene_body_seq)){
						exist.add(gene_body_seq);
					}
					gene_to_gene_body_seq.put(gene_name,exist);
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(gene_body_seq);
					gene_to_gene_body_seq.put(gene_name,exist);
				}

				

			    }//else if(next_type.equals("five_prime_UTR"))


		    }//j
		    
		}//if
	    }//if not ""
	}//w
                    
}
	
 
  String reverse_complement(String a) throws IOException {

	    a = a.toUpperCase();
	    String b = "";

	    for(int i = a.length() -1; i >= 0; i--){

	      String temp = a.substring(i, i+1);

	      if(temp.equals("A")){

		b = b+ "T";

	      }else if(temp.equals("G")){

		b = b+ "C";


	      }else if(temp.equals("C")){

		b = b+ "G";


	      }else if(temp.equals("T")){

		b = b+ "A";
	      }else{

		b = b+ "N";
	      }
	      
	    }

	    return(b);


  }//method

}//inner class dni

static class dna_tf_module{


static LinkedHashMap<String,String> atxg_to_gene_name;
static LinkedHashMap<String,String> gene_name_to_atxg;
 
void  map_atxg_to_gene_name()  throws IOException{

	//System.out.println("map_atxg_to_gene_name():" );

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern colon_pattern = Pattern.compile(":");
	
	atxg_to_gene_name= new LinkedHashMap();
	gene_name_to_atxg=new LinkedHashMap();
 
 	String path = "data/araport.gene.geneName.sort.uniq";
 	
	List<String> lines  = FileUtils.readLines(new File(path));
	
					
	for(int k = 0; k < lines.size(); k++){
				
		String one_line = lines.get(k).toUpperCase();
		String[] space_split = colon_pattern.split(one_line);
		String atxg = space_split[0].trim();
		String gene = space_split[1].trim();
		
		atxg_to_gene_name.put(atxg,gene);
		gene_name_to_atxg.put(gene,atxg);
		
	}
	
}


static LinkedHashMap<String,String> atxg_to_tf_name;
static LinkedHashMap<String,String> tf_name_to_atxg;
 
void  map_atxg_to_tf_name()  throws IOException{

	//System.out.println("map_atxg_to_tf_name():" );

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern comma_pattern = Pattern.compile(",");
	
	atxg_to_tf_name= new LinkedHashMap();
	tf_name_to_atxg=new LinkedHashMap();
 
 	String path = "data/cis_bp_ara_map_tf_gene.csv";
 	//237839
	List<String> lines  = FileUtils.readLines(new File(path));
	//ENSG00000222623:RNU6-1100P
	//NSG00000279928:DDX11L17

					
	for(int k = 0; k < lines.size(); k++){
				
		String one_line = lines.get(k).toUpperCase();
		String[] comma_split = comma_pattern.split(one_line);
		String atxg = comma_split[3].trim();
		String tf = comma_split[1].trim();
		
		atxg_to_gene_name.put(atxg,tf);
		gene_name_to_atxg.put(tf,atxg);
		
	}
	
}



static String[] chr = {"1","2","3","4","5"};

static String[] tf_fam_per_chr ={"ABI3VP1","AP2EREBP","ARF","ARF_ecoli_MP","ARID","BBRBPC","BES1","bHLH","bZIP","BZR","C2C2COlike","C2C2dof","C2C2gata","C2C2YABBY","C2H2","C3H","CAMTA","CCAATHAP3","CPP","DBP","E2FDP","EIL","FAR1","FHA","G2like","GeBP","GRF","HB","HMG","Homeobox","Homeobox_ecoli_HAT2","HSF","LIM","LOBAS2","MADS","mTERF","MYB","MYBrelated","NAC","ND","Orphan","PLATZ","RAV","REM","REMB3","RWPRK","S1Falike","SBP","SRS","TCP","Trihelix","WRKY","zfGRF","ZFHD"};

static String[] tf_per_chr = 
{"AT5G18090","AT5G25475","AT5G60130","FUS3","NGA4","REM16","VRN1","ABR1","AIL7","AT1G01250","AT1G12630","AT1G19210","AT1G22810","AT1G28160","AT1G36060","AT1G44830","AT1G71450","AT1G75490","AT1G77200","AT1G77640","AT2G33710","AT2G44940","AT3G16280","AT3G57600","AT3G60490","AT4G16750","AT4G18450","AT4G28140","AT4G31060","AT4G32800","AT5G18450","AT5G65130","AT5G67000","CBF1","CBF2","CBF3","CBF4","CEJ1","CRF10","CRF4","DDF1","DDF2","DEAR2","DEAR3","DEAR5","DREB19","DREB2","DREB26","ERF1","ERF10","ERF104","ERF105","ERF11","ERF115","ERF13","ERF15","ERF2","ERF3","ERF38","ERF4","ERF48","ERF5","ERF6","ERF7","ERF73","ERF8","ERF9","ESE1","ESE3","LEP","PLT1","PLT3","PUCHI","RAP2.1","Rap2.10","RAP2.11","RAP2.12","RAP2.6","RRTF1","SHN3","TINY","ARF_ecoli_MP","ARF16","ARF2","AT1G04880","AT1G20910","AT1G76110","AT2G17410","AT3G13350","BPC1","BPC5","BPC6","BAM8","bHLH10","bHLH104","bHLH122","bHLH130","bHLH157","bHLH18","bHLH28","bHLH31","bHLH34","bHLH69","bHLH74","bHLH80","BIM1","BIM2","BIM3","PIF7","ABF2","ABI5","AREB3","bZIP16","bZIP18","bZIP28","bZIP3","bZIP42","bZIP43","bZIP44","bZIP48","bZIP50","bZIP52","bZIP53","bZIP68","bZIP69","GBF3","GBF5","GBF6","HY5","TGA1","TGA10","TGA2","TGA3","TGA4","TGA5","TGA6","TGA9","VIP1","AT1G78700","AT4G18890","AT4G36780","BZR1","AT4G27900","AT5G59990","ADOF1","AT1G47655","AT1G64620","AT1G69570","AT2G28810","AT3G45610","AT3G52440","AT4G38000","AT5G02460","AT5G62940","AT5G66940","CDF3","COG1","DAG2","DOF2.4","DOF4.2","DOF4.3","DOF4.5","OBP1","OBP3","OBP4","GATA1","GATA11","GATA12","GATA15","GATA16","GATA19","GATA20","GATA4","ZIM","ZML1","CRC","AT1G14580","AT2G15740","AT2G48100","AT3G46070","AT3G49930","AT3G60580","AT4G26030","AT5G04390","AT5G22890","AT5G22990","AT5G66730","AtIDD11","AZF1","IDD2","IDD4","IDD5","IDD7","JGL","JKD","MGP","NUC","SGR5","STOP1","STZ","TF3A","WIP5","AT1G74370","AT3G12130","AT5G08750","AT5G63260","CDM1","EMB1789","TZF9","U2AF35B","CAMTA1","CAMTA5","HAP3","NFYB4","AT2G20110","SOL1","TCX2","AT3G51470","DEL1","DEL2","E2FA","E2FC","EIL3","EIN3","FAR1","FHA2","AT1G13300","AT1G25550","AT1G49560","AT1G68670","AT2G01060","AT2G03500","AT2G20400","AT2G38300","AT2G40260","AT3G04030","AT3G12730","AT3G24120","AT4G37180","AT5G29000","AT5G45580","KAN2","AT1G66420","AT4G00250","AtGRF6","GRF9","ANL2","ATHB15","ATHB21","ATHB40","ATHB5","ATHB53","HDG7","LMI1","PHV","WOX11","3XHMGBOX1","Homeobox_ecoli_HAT2","ATHB13","ATHB18","ATHB20","ATHB6","ATHB7","EDT1","HAT1","HAT2","HAT5","HDG1","PDF2","WUS1","HSF21","HSF3","HSF6","HSF7","HSFA1E","HSFA6A","HSFA6B","HSFB3","HSFB4","HSFC1","WLIM2A","AS2","ASL18","LBD13","LBD18","LBD19","LBD2","LBD23","LOB","AGL13","AGL15","AGL16","AGL25","AGL42","AGL55","AGL6","AGL63","FEM111","SVP","AT5G23930","AT1G18960","AT1G19000","AT1G49010","AT1G72740","AT1G74840","AT2G38090","AT3G09600","AT3G10113","AT3G10580","AT3G11280","AT4G01280","AT4G12670","AT5G05790","AT5G08520","AT5G47390","AT5G52660","AT5G56840","AT5G58900","AT5G61620","EPR1","LCL1","LHY1","RVE1","TBP3","TRP1","TRP2","ATY13","ATY19","BOS1","MS188","MYB1","MYB10","MYB101","MYB105","MYB107","MYB113","MYB116","MYB118","MYB119","MYB121","MYB13","MYB17","MYB23","MYB27","MYB30","MYB33","MYB39","MYB3R1","MYB3R4","MYB3R5","MYB4","MYB40","MYB41","MYB43","MYB44","MYB49","MYB51","MYB52","MYB55","MYB56","MYB57","MYB58","MYB60","MYB61","MYB62","MYB63","MYB65","MYB67","MYB70","MYB73","MYB74","MYB77","MYB81","MYB83","MYB88","MYB92","MYB93","MYB94","MYB96","MYB98","MYB99","ANAC004","ANAC005","ANAC011","ANAC013","ANAC016","ANAC017","ANAC020","ANAC028","ANAC034","ANAC038","ANAC042","ANAC045","ANAC046","ANAC047","ANAC050","ANAC053","ANAC055","ANAC057","ANAC058","ANAC062","ANAC070","ANAC071","ANAC075","ANAC079","ANAC083","ANAC087","ANAC092","ANAC094","ANAC096","ANAC103","AT1G19040","AT3G12910","ATAF1","CUC1","CUC2","CUC3","NAC2","NAM","NAP","NST1","NTL8","NTM1","NTM2","SMB","SND2","SND3","VND1","VND2","VND3","VND4","VND6","AGL95","AT1G63040","AT2G28920","FRS9","AT1G23810","AT1G24250","BBX31","AT2G01818","RAV1","AT2G31460","REM19","NLP7","RKD2","AT3G09735","SPL1","SPL11","SPL13","SPL14","SPL15","SPL3","SPL5","SPL9","SRS7","AT1G69690","AT1G72010","AT2G45680","AT5G08330","PTF1","TCP1","TCP16","TCP17","TCP20","TCP24","TCP3","TCP7","AT1G76870","AT1G76880","AT2G33550","AT3G10030","AT3G14180","AT3G25990","AT3G58630","AT5G05550","AT5G47660","GT1","GT2","GT3a","GTL1","WRKY11","WRKY14","WRKY15","WRKY17","WRKY18","WRKY20","WRKY21","WRKY22","WRKY24","WRKY25","WRKY26","WRKY27","WRKY28","WRKY29","WRKY3","WRKY30","WRKY31","WRKY33","WRKY40","WRKY42","WRKY43","WRKY45","WRKY46","WRKY47","WRKY50","WRKY55","WRKY59","WRKY6","WRKY65","WRKY7","WRKY70","WRKY71","WRKY75","WRKY8","AT3G42860","ATHB23","ATHB24","ATHB25","AtHB32","ATHB33","ATHB34"};

static LinkedHashMap<String,String> tf_to_tf_fam;
static LinkedList<String>[] chip_tf_per_chr;   
static String[][] chip_tf_per_chr_str;


void link_tf_to_fam() throws IOException{
/*
ABI3VP1_tnt_AT5G18090
ABI3VP1_tnt_AT5G25475
ABI3VP1_tnt_AT5G60130
*/
	System.out.println("link_tf_to_fam:" );

	Pattern underbar_pattern = Pattern.compile("_");

	tf_to_tf_fam = new LinkedHashMap();
	
	chip_tf_per_chr = new LinkedList[tf_fam_per_chr.length];
	chip_tf_per_chr_str = new String[tf_fam_per_chr.length][];

	for(int j = 0; j < tf_fam_per_chr.length; j++){

		chip_tf_per_chr[j] = new LinkedList();

	}
	


	String dir_path = "big/chip_analysis/result_csv/tf_sort_uniq/";

	String file_path = dir_path + "comb.chr1.tf1.sort.uniq";
	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 0; q < lines.size(); q++){
		
		String one_line = lines.get(q).trim();
		String[] bar_split = underbar_pattern.split(one_line);
			
		String fam = bar_split[0].trim();
		String tf = bar_split[2].trim();
		tf_to_tf_fam.put(tf,fam);
			
		for(int j = 0; j < tf_fam_per_chr.length; j++){

			if(fam.equals(tf_fam_per_chr[j])){

				chip_tf_per_chr[j].add(tf);
				break;
			}

		}

			
	}//q
	
	for(int q = 0; q < tf_fam_per_chr.length; q++){

		chip_tf_per_chr_str[q] = new String[chip_tf_per_chr[q].size()];

		for(int w = 0; w < chip_tf_per_chr[q].size();w++){

			chip_tf_per_chr_str[q][w] = chip_tf_per_chr[q].get(w);
		}
	}
		
	

}//method

    
//retrieve tf_comb chip-seq score

static LinkedList<Double>[][][] retrieve_tf_comb_chip_seq_score(String file_path) throws IOException{

	System.out.println("retrieve_tf_comb_chip_seq_score():" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][] tf_score = new LinkedList[tf_per_chr.length][tf_per_chr.length][5];
	
	
	for(int j = 0; j < tf_per_chr.length; j++){

		for(int k = 0; k < tf_per_chr.length; k++){

			for(int l = 0; l < 5; l++){
					tf_score[j][k][l] = new LinkedList();
			}

		}
	}
	

	
	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
		/*
		ARF_tnt_ARF16-Trihelix_tnt_AT3G25990,,ARF_tnt_ARF16,Trihelix_tnt_AT3G25990,25,55,367,0.17596473175581934,91.8563386845009
Trihelix_tnt_AT3G25990-ARF_tnt_ARF16,,Trihelix_tnt_AT3G25990,ARF_tnt_ARF16,25,367,55,0.17596473175581934,91.8563386845009

		*/
		String one_line = lines.get(q).trim();

		//System.out.println(one_line);
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1_split = comma_split[2].trim();
		String tf2_split = comma_split[3].trim();

		String[] bar_split1 = underbar_pattern.split(tf1_split);
			
		String tf1 = bar_split1[2].trim();

		String[] bar_split2 = underbar_pattern.split(tf2_split);
			
		String tf2 = bar_split2[2].trim();
		
		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}

		boolean tf_proceed = false;

		if(tf_int1 != -1 && tf_int2 != -1){

			tf_proceed = true;

		}

		if(tf_proceed){

			Double temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			tf_score[tf_int1][tf_int2][0].add(temp_double);

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			tf_score[tf_int1][tf_int2][1].add(temp_double);

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			tf_score[tf_int1][tf_int2][2].add(temp_double);

			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			tf_score[tf_int1][tf_int2][3].add(temp_double);

			temp_double = new Double(0);
			if(comma_split[8].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[8].trim());
			}
			tf_score[tf_int1][tf_int2][4].add(temp_double);

		}

			
	}//q

	return(tf_score);	
	
}//method

   
//retrieve tf_comb chip-seq score

static LinkedList<Double>[][][] retrieve_tf_comb_chip_seq_score_fam(String file_path) throws IOException{

	System.out.println("retrieve_tf_comb_chip_seq_score_fam():" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][] tf_fam_score = new LinkedList[tf_fam_per_chr.length][tf_fam_per_chr.length][5];

	for(int j = 0; j < tf_fam_per_chr.length; j++){

		for(int k = 0; k < chip_tf_per_chr_str.length; k++){

			for(int l = 0; l < 5; l++){

				tf_fam_score[j][k][l] = new LinkedList();
			}
			
		}
	}

	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
		/*
		ARF_tnt_ARF16-Trihelix_tnt_AT3G25990,,ARF_tnt_ARF16,Trihelix_tnt_AT3G25990,25,55,367,0.17596473175581934,91.8563386845009
Trihelix_tnt_AT3G25990-ARF_tnt_ARF16,,Trihelix_tnt_AT3G25990,ARF_tnt_ARF16,25,367,55,0.17596473175581934,91.8563386845009

		*/
		String one_line = lines.get(q).trim();

		//System.out.println(one_line);
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1_split = comma_split[2].trim();
		String tf2_split = comma_split[3].trim();

		String[] bar_split1 = underbar_pattern.split(tf1_split);
			
		String fam1 = bar_split1[0].trim();
		String tf1 = bar_split1[2].trim();

		String[] bar_split2 = underbar_pattern.split(tf2_split);
			
		String fam2 = bar_split2[0].trim();
		String tf2 = bar_split2[2].trim();
			
		int fam_int1 = -1;
		for(int j = 0; j < tf_fam_per_chr.length; j++){

			if(fam1.equals(tf_fam_per_chr[j])){

				fam_int1 =j;
				break;
			}

		}

		int fam_int2 = -1;
		for(int j = 0; j < tf_fam_per_chr.length; j++){

			if(fam2.equals(tf_fam_per_chr[j])){

				fam_int2 =j;
				break;
			}

		}

		boolean fam_proceed = false;

		if(fam_int1 != -1 && fam_int2 != -1){

			fam_proceed = true;

		}

		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}

		boolean tf_proceed = false;

		if(tf_int1 != -1 && tf_int2 != -1){

			tf_proceed = true;

		}

		if(tf_proceed){

			Double temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			
			if(fam_proceed){
				tf_fam_score[fam_int1][fam_int2][0].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			
			if(fam_proceed){
				tf_fam_score[fam_int1][fam_int2][1].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			
			if(fam_proceed){
				tf_fam_score[fam_int1][fam_int2][2].add(temp_double);
			}


			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			
			if(fam_proceed){
				tf_fam_score[fam_int1][fam_int2][3].add(temp_double);
			}


			temp_double = new Double(0);
			if(comma_split[8].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[8].trim());
			}
			
			if(fam_proceed){
				tf_fam_score[fam_int1][fam_int2][4].add(temp_double);
			}
		}

			
	}//q

	return(tf_fam_score);	
	
}//method


//if tf has not match, use this method
//input non-match tf
//returns match-tf
static String get_match_tf(String tf){

      String new_tf = "";
      if(tf.matches("^[A-Z]+[0]+[0-9]+.*")){
		new_tf = "";
		int letter_part = -1;
		
		for(int n = 0; n < tf.length(); n++){

			String cur = tf.substring(n,n+1);
			if(cur.equals("0")){

				new_tf = tf.substring(0,n);
				letter_part = n;
				break;
			}
		}
		
		int digit_part = -1;
		for(int n = letter_part; n < tf.length(); n++){

			String cur = tf.substring(n,n+1);
			if(!cur.equals("0")){

				new_tf = new_tf + tf.substring(n, tf.length());
				digit_part = n;
				
				break;
			}
		}
	}

	if(!new_tf.equals("")){
		tf = new_tf;
	}

	boolean check = false;
	for(int i = 0; i < tf_per_chr.length; i++){

		if(tf.equals(tf_per_chr[i])){

			check = true;
			break;
		}
	}


	if(!check){

		if(tf.matches("^[A-Z]+[0-9]+[A-Z]+.*")){

			for(int j = tf.length()-1; j >= 0; j--){

				String cur = tf.substring(j,j+1);

				if(cur.matches("[0-9]")){

					tf = tf.substring(0,j+1);
					break;
				}
			}
		}


	}
	return(tf);
}

//atac analysis. orientation

static String[] sign_comb = {"++","+-","-+","--"};

//retrieve atac orientation score

LinkedList<Double>[][][][] retrieve_atac_orien_score(String file_path) throws IOException{

	System.out.println("retrieve_atac_orien_score():" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][][] atac_orien_tf_score = new LinkedList[tf_per_chr.length][][][];

	for(int j = 0; j < tf_per_chr.length; j++){

		atac_orien_tf_score[j] = new LinkedList[tf_per_chr.length][][];

		for(int k = 0; k < tf_per_chr.length; k++){

			atac_orien_tf_score[j][k] = new LinkedList[sign_comb.length][];

			for(int l = 0; l < sign_comb.length; l++){

				atac_orien_tf_score[j][k][l] = new LinkedList[5];

				for(int m = 0; m < 5; m++){
					atac_orien_tf_score[j][k][l][m] = new LinkedList();
				}
				

			}

		}
	}
	

	
	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
	,TF1,TF2,TF1_TF2_count,TF1_count,TF2_count,cosine,zscore
DREB2D(+)-ERF057(+),DREB2D(+),ERF057(+),114,125,128,0.9012491331479882,151.2698629224322

	*/
		String one_line = lines.get(q).trim();
		//System.out.println("849:" + one_line);
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1_split = comma_split[1].trim();
		String tf2_split = comma_split[2].trim();

		String sign1 = "";
		if(tf1_split.contains("+")){

			sign1 = "+";

		}else if(tf1_split.contains("-")){
			sign1 = "-";

		}
		String tf1 = tf1_split.substring(0,tf1_split.length()-3);

		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;

		String sign2 = "";
		if(tf2_split.contains("+")){

			sign2 = "+";

		}else if(tf2_split.contains("-")){
			sign2 = "-";

		}
		String tf2 = tf2_split.substring(0,tf2_split.length()-3);

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;
			
		String cur_sign = sign1+sign2;
		int sign_int = -1;
		for(int j = 0; j < sign_comb.length; j++){

			if(cur_sign.equals(sign_comb[j])){

				sign_int =j;
				break;
			}

		}

		System.out.println("atac tf:" +tf1 + ":" + tf2);
		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
		System.out.println("atac tf_int1:" +tf_int1+ ":" + tf_int2);
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		
 		if(tf_proceed){
			Double temp_double = new Double(0);
			if(comma_split[3].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[3].trim());
			}
			atac_orien_tf_score[tf_int1][tf_int2][sign_int][0].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			atac_orien_tf_score[tf_int1][tf_int2][sign_int][1].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			atac_orien_tf_score[tf_int1][tf_int2][sign_int][2].add(temp_double);
				
			

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			atac_orien_tf_score[tf_int1][tf_int2][sign_int][3].add(temp_double);
				
			

			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			atac_orien_tf_score[tf_int1][tf_int2][sign_int][4].add(temp_double);

			
		}

			
	}//q

	return(atac_orien_tf_score);	
	
}//method

//retrieve atac orientation score

static LinkedList<Double>[][][][] retrieve_atac_orien_score_fam(String file_path) throws IOException{

	System.out.println("retrieve_atac_orien_score():" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][][] atac_orien_tf_fam_score = new LinkedList[chip_tf_per_chr_str.length][][][];

	for(int j = 0; j < chip_tf_per_chr_str.length; j++){

		atac_orien_tf_fam_score[j] = new LinkedList[chip_tf_per_chr_str.length][][];

		for(int k = 0; k < chip_tf_per_chr_str.length; k++){

			atac_orien_tf_fam_score[j][k] = new LinkedList[sign_comb.length][];

			for(int l = 0; l < sign_comb.length; l++){
	
				atac_orien_tf_fam_score[j][k][l] = new LinkedList[5];
			
				for(int m = 0; m < 5; m++){


					atac_orien_tf_fam_score[j][k][l][m] = new LinkedList();

				}
			}
		}
	}
	

	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
	,TF1,TF2,TF1_TF2_count,TF1_count,TF2_count,cosine,zscore
DREB2D(+)-ERF057(+),DREB2D(+),ERF057(+),114,125,128,0.9012491331479882,151.2698629224322

	*/
		String one_line = lines.get(q).trim();
		//System.out.println("849:" + one_line);
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1_split = comma_split[1].trim();
		String tf2_split = comma_split[2].trim();

		String sign1 = "";
		if(tf1_split.contains("+")){

			sign1 = "+";

		}else if(tf1_split.contains("-")){
			sign1 = "-";

		}
		String tf1 = tf1_split.substring(0,tf1_split.length()-3);

		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;

		String sign2 = "";
		if(tf2_split.contains("+")){

			sign2 = "+";

		}else if(tf2_split.contains("-")){
			sign2 = "-";

		}
		String tf2 = tf2_split.substring(0,tf2_split.length()-3);

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;
			
		String cur_sign = sign1+sign2;
		int sign_int = -1;
		for(int j = 0; j < sign_comb.length; j++){

			if(cur_sign.equals(sign_comb[j])){

				sign_int =j;
				break;
			}

		}

		System.out.println("atac tf:" +tf1 + ":" + tf2);
		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
		System.out.println("atac tf_int1:" +tf_int1+ ":" + tf_int2);
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		String fam1 = "";
		String fam2 = "";
		boolean fam_proceed = false;

		if(tf_proceed){
			
			fam1 = tf_to_tf_fam.get(tf1);
			fam2 = tf_to_tf_fam.get(tf2);

			if(fam1 != null && fam2 != null){

				fam_proceed = true;
			}
				
		}

		int fam_int1 = -1;
		int fam_int2 = -1;

		if(fam_proceed){
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam1.equals(tf_fam_per_chr[j])){

					fam_int1 =j;
					break;
				}

			}

				
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam2.equals(tf_fam_per_chr[j])){

					fam_int2 =j;
					break;
				}

			}
		}
			

 		System.out.println("atac fam_int1:" +fam_int1+ ":" +fam_int2);
		if(tf_proceed){
			Double temp_double = new Double(0);
			if(comma_split[3].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[3].trim());
			}
			if(fam_proceed){
				atac_orien_tf_fam_score[fam_int1][fam_int2][sign_int][0].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			if(fam_proceed){
				atac_orien_tf_fam_score[fam_int1][fam_int2][sign_int][1].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			if(fam_proceed){
				atac_orien_tf_fam_score[fam_int1][fam_int2][sign_int][2].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			if(fam_proceed){
				atac_orien_tf_fam_score[fam_int1][fam_int2][sign_int][3].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			if(fam_proceed){
				atac_orien_tf_fam_score[fam_int1][fam_int2][sign_int][4].add(temp_double);	
			}
		}

			
	}//q

	return(atac_orien_tf_fam_score);	
	
}//method

    
//retrieve tf direction atac score

static LinkedList<Double>[][][] retrieve_tf_direction_score(String file_path) throws IOException{

	System.out.println("retrieve_tf_direction_score:" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][] tf_direction_score = new LinkedList[tf_per_chr.length][tf_per_chr.length][7];
	
	tf_direction_score= new LinkedList[tf_per_chr.length][][];

	for(int j = 0; j < tf_per_chr.length; j++){

		for(int k = 0; k < tf_per_chr.length; k++){

			for(int l = 0; l < 7; l++){
				tf_direction_score[j][k][l] = new LinkedList();
			}

		}

	}
	


	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
		,TF1,TF2,TF1_TF2_count,TF1-TF2,TF2-TF1,convergent,divergent,std,pvalue
1-1,1,1,6728,0.3329369797859691,0.3329369797859691,0.16810344827586207,0.16602259215219975,0.09577114291402651,3.431637383616119e-160
MNB1A-MNB1A,MNB1A,MNB1A,3408,0.35944835680751175,0.35944835680751175,0.14847417840375587,0.13262910798122066,0.12654551797647776,1.2623060834837178e-141
PBF-PBF,PBF,PBF,3408,0.35944835680751175,0.35944835680751175,0.14847417840375587,0.13262910798122066,0.12654551797647776,1.2623060834837178e-141

	*/
		String one_line = lines.get(q).trim();
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1= comma_split[1].trim();
		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;
		String tf2 = comma_split[2].trim();

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;
		System.out.println("direction tf:" +tf1+ ":" +tf2);
		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
		System.out.println("direction tf_int1:" +tf_int1+ ":" +tf_int2);
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		
		if(tf_proceed){

			Double temp_double = new Double(0);
			if(comma_split[3].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[3].trim());
			}
			tf_direction_score[tf_int1][tf_int2][0].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			tf_direction_score[tf_int1][tf_int2][1].add(temp_double);


			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			tf_direction_score[tf_int1][tf_int2][2].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			tf_direction_score[tf_int1][tf_int2][3].add(temp_double);

			
			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			tf_direction_score[tf_int1][tf_int2][4].add(temp_double);

			
			temp_double = new Double(0);
			if(comma_split[8].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[8].trim());
			}
			tf_direction_score[tf_int1][tf_int2][5].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[9].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[9].trim());
			}
			tf_direction_score[tf_int1][tf_int2][6].add(temp_double);

			
		}

			
	}//q

	return(tf_direction_score);	
	
}//method

//retrieve tf direction atac score

static LinkedList<Double>[][][] retrieve_tf_direction_score_fam(String file_path) throws IOException{

	System.out.println("retrieve_tf_direction_score_fam:" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	
	LinkedList<Double>[][][] tf_direction_fam_score = new LinkedList[chip_tf_per_chr_str.length][chip_tf_per_chr_str.length][7];

	for(int j = 0; j < chip_tf_per_chr_str.length; j++){

		for(int k = 0; k < chip_tf_per_chr_str.length; k++){

			for(int l = 0; l < 7; l++){

				tf_direction_fam_score[j][k][l] = new LinkedList();
			}
		}
	}
	

	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
		,TF1,TF2,TF1_TF2_count,TF1-TF2,TF2-TF1,convergent,divergent,std,pvalue
1-1,1,1,6728,0.3329369797859691,0.3329369797859691,0.16810344827586207,0.16602259215219975,0.09577114291402651,3.431637383616119e-160
MNB1A-MNB1A,MNB1A,MNB1A,3408,0.35944835680751175,0.35944835680751175,0.14847417840375587,0.13262910798122066,0.12654551797647776,1.2623060834837178e-141
PBF-PBF,PBF,PBF,3408,0.35944835680751175,0.35944835680751175,0.14847417840375587,0.13262910798122066,0.12654551797647776,1.2623060834837178e-141

	*/
		String one_line = lines.get(q).trim();
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1= comma_split[1].trim();
		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;
		String tf2 = comma_split[2].trim();

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;
		System.out.println("direction tf:" +tf1+ ":" +tf2);
		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
		System.out.println("direction tf_int1:" +tf_int1+ ":" +tf_int2);
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		String fam1 = "";
		String fam2 = "";
		boolean fam_proceed = false;

		if(tf_proceed){
		
			fam1 = tf_to_tf_fam.get(tf1);
			fam2 = tf_to_tf_fam.get(tf2);

			if(fam1 != null && fam2 != null){

				fam_proceed = true;
			}
				
		}
		System.out.println("direction fam1:" +fam1+ ":" +fam2);
		int fam_int1 = -1;
		int fam_int2 = -1;

		if(fam_proceed){
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam1.equals(tf_fam_per_chr[j])){

					fam_int1 =j;
					break;
				}

			}

				
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam2.equals(tf_fam_per_chr[j])){

					fam_int2 =j;
					break;
				}

			}

			if(!(fam_int1 != -1 && fam_int2 != -1)){

				fam_proceed = false;
			}
		}

		System.out.println("direction fam_int1:" +fam_int1+ ":" +fam_int2);
		if(tf_proceed){

			Double temp_double = new Double(0);
			if(comma_split[3].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[3].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][0].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][0].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][1].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][2].add(temp_double);
			}


			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][3].add(temp_double);
			}


			temp_double = new Double(0);
			if(comma_split[8].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[8].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][4].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[9].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[9].trim());
			}
			if(fam_proceed){
				tf_direction_fam_score[fam_int1][fam_int2][4].add(temp_double);
			}
		}

			
	}//q

	return(tf_direction_fam_score);	
	
}//method


//retrieve tf same_opposite atac score

static LinkedList<Double>[][][] retrieve_tf_same_opposite_score(String file_path) throws IOException{

	System.out.println("retrieve_tf_same_opposite_score:" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][] tf_same_opposite_score = new LinkedList[tf_per_chr.length][tf_per_chr.length][5];
	
	for(int j = 0; j < tf_per_chr.length; j++){

		for(int k = 0; k < tf_per_chr.length; k++){

			for(int l = 0; l < 5; l++){
				tf_same_opposite_score[j][k][l] = new LinkedList();
			}

		}

	}


	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
		,TF1,TF2,TF1_TF2_count,same,opposite,std,pvalue
DREB2D-ERF057,DREB2D,ERF057,129,0.9844961240310077,0.015503875968992248,0.6851809895218484,3.5897266476400753e-28
DREB2D-ERF014,DREB2D,ERF014,113,1.0,0.0,0.7071067811865476,2.157748546047867e-26
ERF014-ERF057,ERF014,ERF057,115,0.991304347826087,0.008695652173913044,0.6948092719485206,5.814158857566564e-26
	*/
		String one_line = lines.get(q).trim();
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1= comma_split[1].trim();
		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;
		String tf2 = comma_split[2].trim();

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;

		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
			
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		
		if(tf_proceed){

			Double temp_double = new Double(0);
			if(comma_split[3].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[3].trim());
			}
			tf_same_opposite_score[tf_int1][tf_int2][0].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			tf_same_opposite_score[tf_int1][tf_int2][1].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			tf_same_opposite_score[tf_int1][tf_int2][2].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			tf_same_opposite_score[tf_int1][tf_int2][3].add(temp_double);

			

			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			tf_same_opposite_score[tf_int1][tf_int2][4].add(temp_double);

			
				
		}

			
	}//q
	
	return(tf_same_opposite_score);
}//method

//retrieve tf_same_opposite_fam_score

static LinkedList<Double>[][][] retrieve_tf_same_opposite_score_fam(String file_path) throws IOException{

	System.out.println("retrieve_tf_same_opposite_score_fam:" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][] tf_same_opposite_fam_score = new LinkedList[chip_tf_per_chr_str.length][][];

	for(int j = 0; j < chip_tf_per_chr_str.length; j++){

		tf_same_opposite_fam_score[j] = new LinkedList[chip_tf_per_chr_str.length][];

		for(int k = 0; k < chip_tf_per_chr_str.length; k++){

			tf_same_opposite_fam_score[j][k] = new LinkedList[5];

			for(int l = 0; l < 5; l++){

				tf_same_opposite_fam_score[j][k][l] = new LinkedList();
			}
		}
	}
	

	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
		,TF1,TF2,TF1_TF2_count,same,opposite,std,pvalue
DREB2D-ERF057,DREB2D,ERF057,129,0.9844961240310077,0.015503875968992248,0.6851809895218484,3.5897266476400753e-28
DREB2D-ERF014,DREB2D,ERF014,113,1.0,0.0,0.7071067811865476,2.157748546047867e-26
ERF014-ERF057,ERF014,ERF057,115,0.991304347826087,0.008695652173913044,0.6948092719485206,5.814158857566564e-26
	*/
		String one_line = lines.get(q).trim();
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1= comma_split[1].trim();
		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;
		String tf2 = comma_split[2].trim();

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;

		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
			
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		String fam1 = "";
		String fam2 = "";
		boolean fam_proceed = false;

		if(tf_proceed){
			
			fam1 = tf_to_tf_fam.get(tf1);
			fam2 = tf_to_tf_fam.get(tf2);

			if(fam1 != null && fam2 != null){

				fam_proceed = true;
			}
				
		}

		int fam_int1 = -1;
		int fam_int2 = -1;

		if(fam_proceed){
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam1.equals(tf_fam_per_chr[j])){

					fam_int1 =j;
					break;
				}

			}

				
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam2.equals(tf_fam_per_chr[j])){

					fam_int2 =j;
					break;
				}

			}
		}

			
		if(tf_proceed){

			Double temp_double = new Double(0);
			if(comma_split[3].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[3].trim());
			}
			
			if(fam_proceed){
				tf_same_opposite_fam_score[fam_int1][fam_int2][0].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[4].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[4].trim());
			}
			
			if(fam_proceed){
				tf_same_opposite_fam_score[fam_int1][fam_int2][0].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[5].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[5].trim());
			}
			
			if(fam_proceed){
				tf_same_opposite_fam_score[fam_int1][fam_int2][1].add(temp_double);
			}

			temp_double = new Double(0);
			if(comma_split[6].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[6].trim());
			}
			
			if(fam_proceed){
				tf_same_opposite_fam_score[fam_int1][fam_int2][2].add(temp_double);
			}


			temp_double = new Double(0);
			if(comma_split[7].trim().equals("inf")){

				temp_double = Double.POSITIVE_INFINITY;
			}else{
				temp_double = new Double(comma_split[7].trim());
			}
			
			if(fam_proceed){
				tf_same_opposite_fam_score[fam_int1][fam_int2][3].add(temp_double);
			}


				
		}

			
	}//q
	
	return(tf_same_opposite_fam_score);
}//method


//retrieve tf distance atac score

static LinkedList<Double>[][][] retrieve_tf_distance_score(String file_path) throws IOException{

	System.out.println("retrieve_tf_distance_score:" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	LinkedList<Double>[][][] tf_distance_score = new LinkedList[tf_per_chr.length][tf_per_chr.length][10];

	for(int j = 0; j < tf_per_chr.length; j++){

		for(int k = 0; k < tf_per_chr.length; k++){

			for(int l = 0; l < 101; l++){
				tf_distance_score[j][k][l] = new LinkedList();
			}

		}

	}
	


	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
		,TF1,TF2,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100

PBF-MNB1A,PBF,MNB1A,62,100,68,48,56,56,54,46,38,52,50,70,36,70,42,40,72,50,28,66,46,36,46,42,48,54,48,46,48,60,52,56,76,40,42,52,40,26,54,26,42,36,28,28,46,62,26,38,52,60,38,46,46,50,28,32,52,40,38,44,32,46,50,30,24,48,44,56,28,32,32,24,20,26,36,44,34,30,46,38,40,42,26,26,48,30,38,34,52,38,32,30,44,66,34,36,20,40,26,26,44
	*/
		String one_line = lines.get(q).trim();
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1= comma_split[1].trim();
		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;
		String tf2 = comma_split[2].trim();

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;

		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
			
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

			
		if(tf_proceed){
	
			for(int z = 3; z < comma_split.length; z++){

				Double temp_double = new Double(0);
				if(comma_split[z].trim().equals("inf")){

					temp_double = Double.POSITIVE_INFINITY;
				}else{
					temp_double = new Double(comma_split[z].trim());
				}
				tf_distance_score[tf_int1][tf_int2][z-3].add(temp_double);

				
			}

				

				
		}

			
	}//q
	
	return(tf_distance_score);
}//method

//retrieve tf distance atac score

static LinkedList<Double>[][][] retrieve_tf_distance_score_fam(String file_path) throws IOException{

	System.out.println("retrieve_tf_distance_score:" );
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern underbar_pattern = Pattern.compile("_");

	
	LinkedList<Double>[][][] tf_distance_fam_score = new LinkedList[chip_tf_per_chr_str.length][][];

	for(int j = 0; j < chip_tf_per_chr_str.length; j++){

		tf_distance_fam_score[j] = new LinkedList[chip_tf_per_chr_str.length][];

		for(int k = 0; k < chip_tf_per_chr_str.length; k++){

			tf_distance_fam_score[j][k] = new LinkedList[101];

			for(int l = 0; l < 101; l++){

				tf_distance_fam_score[j][k][l] = new LinkedList();
			}
		}
	}
	

	List<String> lines  = FileUtils.readLines(new File(file_path));
	
	for(int q= 1; q < lines.size(); q++){
	/*
		,TF1,TF2,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100

PBF-MNB1A,PBF,MNB1A,62,100,68,48,56,56,54,46,38,52,50,70,36,70,42,40,72,50,28,66,46,36,46,42,48,54,48,46,48,60,52,56,76,40,42,52,40,26,54,26,42,36,28,28,46,62,26,38,52,60,38,46,46,50,28,32,52,40,38,44,32,46,50,30,24,48,44,56,28,32,32,24,20,26,36,44,34,30,46,38,40,42,26,26,48,30,38,34,52,38,32,30,44,66,34,36,20,40,26,26,44
	*/
		String one_line = lines.get(q).trim();
		String[] comma_split = comma_pattern.split(one_line);
			
		String tf1= comma_split[1].trim();
		boolean check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf1.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		String new_tf = tf1;
		if(!check){

			new_tf = get_match_tf(tf1);
		}
		tf1 = new_tf;
		String tf2 = comma_split[2].trim();

		check = false;
		for(int d = 0; d < tf_per_chr.length; d++){

			if(tf2.equals(tf_per_chr[d])){

				check = true;
				break;
			}
		}

		new_tf = tf2;
		if(!check){

			new_tf = get_match_tf(tf2);
		}
		tf2 = new_tf;

		int tf_int1 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf1.equals(tf_per_chr[j])){

				tf_int1 =j;
				break;
			}

		}

		int tf_int2 = -1;
		for(int j = 0; j < tf_per_chr.length; j++){

			if(tf2.equals(tf_per_chr[j])){

				tf_int2 =j;
				break;
			}

		}
			
		boolean tf_proceed = false;
		if(tf_int1 != -1 && tf_int2 != -1){


			tf_proceed = true;
		}

		String fam1 = "";
		String fam2 = "";
		boolean fam_proceed = false;

		if(tf_proceed){
			
			fam1 = tf_to_tf_fam.get(tf1);
			fam2 = tf_to_tf_fam.get(tf2);

			if(fam1 != null && fam2 != null){

				fam_proceed = true;
			}
				
		}

		int fam_int1 = -1;
		int fam_int2 = -1;

		if(fam_proceed){
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam1.equals(tf_fam_per_chr[j])){

					fam_int1 =j;
					break;
				}

			}

				
			for(int j = 0; j < tf_fam_per_chr.length; j++){

				if(fam2.equals(tf_fam_per_chr[j])){

					fam_int2 =j;
					break;
				}

			}
		}

			
		if(tf_proceed){
	
			for(int z = 3; z < comma_split.length; z++){

				Double temp_double = new Double(0);
				if(comma_split[z].trim().equals("inf")){

					temp_double = Double.POSITIVE_INFINITY;
				}else{
					temp_double = new Double(comma_split[z].trim());
				}
				if(fam_proceed){
					tf_distance_fam_score[fam_int1][fam_int2][z-3].add(temp_double);
				}
			}

				

				
		}

			
	}//q
	
	return(tf_distance_fam_score);
}//method

}//inner class dna_tf 

static class DNA_struc_metal_binding{

//cite:
LinkedList<String> ara_gquad_ortho;
void get_ara_gquad_ortho() throws IOException{
	
	System.out.println("get_ara_gquad_ortho()" );

	ara_gquad_ortho = new LinkedList();

	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/ara_gene_body_gquad_ortholog";
	
	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
							
		String line = lines.get(q).trim();
		ara_gquad_ortho.add(line);
	}
}


String[] ara_target_with_g3 = {"AT5G63290","AT5G05700","AT5G12180","AT5G15550"};

String[] gq_type = {"G2L1-2","G2L1-4","G2L1-9","G3Bulge","G3L1-15","G3VL1-9"};

LinkedHashMap<String,Boolean[]> gene_to_k_dependent_qquad;

void map_gene_to_k_dependent_qquad() throws IOException {
System.out.println("map_gene_to_k_dependent_qquad():" );

	gene_to_k_dependent_qquad = new LinkedHashMap();
		
/*
AT5G01712.1,141,114,1868,1566,1578,3'UTR,GGTAGGGGGACGG,G2L1-2
AT1G22280.2,221,600,721,1398,1430,3'UTR,GGGGGGTAGTGGGTGGTGTTTTGAGGTAAGGGG,G3Bulge
*/


	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern semi_colon_pattern = Pattern.compile(";");
	Pattern digit_pattern = Pattern.compile("\\d");

	String path = "data/K_plus_dependent_g_quad.csv";
	List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){		

		String line = lines.get(q).trim();

		System.out.println(line);
		String[] tab_split = comma_pattern.split(line);

		if(tab_split.length == 9){
				
			String cur_gene = tab_split[0].trim();
			cur_gene = cur_gene.substring(0,cur_gene.indexOf("."));
			String cur_region = tab_split[6].trim();
			
			if(cur_region.equals("CDS")){
				String cur_type = tab_split[8].trim();
				int type_int = -1;
				
				for(int i = 0; i < gq_type.length; i++){
				
					if(cur_type.equals(gq_type[i])){

						type_int = i;
						break;
					}
				}

				if(type_int == -1 ){
				
				}else{
				
					if(gene_to_k_dependent_qquad.containsKey(cur_gene)){
				
						Boolean[] exist = gene_to_k_dependent_qquad.get(cur_gene);
						exist[type_int] = new Boolean(true);
							
						gene_to_k_dependent_qquad.put(cur_gene,exist);
							
					}else{
						Boolean[] exist = new Boolean[gq_type.length];
						for(int r = 0; r < gq_type.length; r++){
						
							if(r == type_int){
								exist[r] = new Boolean(true);
							}else{

								exist[r] = new Boolean(false);
							}
						}
						gene_to_k_dependent_qquad.put(cur_gene,exist);
					}
				
				}
					
			}//7
		}//if
			
			
	}//q
	
}//method

//get g-quad forming genes (only the genes forming g-quad in cds)
LinkedHashMap<String,Boolean[]> gene_to_k_PDS_dependent_qquad;
void map_gene_to_k_PDS_dependent_qquad() throws IOException {
System.out.println("map_gene_to_k_PDS_dependent_qquad():" );

	gene_to_k_PDS_dependent_qquad = new LinkedHashMap();
		
/*
AT5G01712.1,141,114,1868,1566,1578,3'UTR,GGTAGGGGGACGG,G2L1-2
AT1G22280.2,221,600,721,1398,1430,3'UTR,GGGGGGTAGTGGGTGGTGTTTTGAGGTAAGGGG,G3Bulge
*/


	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern semi_colon_pattern = Pattern.compile(";");
	Pattern digit_pattern = Pattern.compile("\\d");

	String path = "data/K_PDS_dependent_g_quad.csv";
	List<String> lines  = FileUtils.readLines(new File(path));
	/*
	AT1G01100.2,83,339,209,282,292,CDS,GGTGGAGGTGG,G2L1-2
AT1G01720.1,141,870,305,940,951,CDS,GGGAGGAGGAGG,G2L1-2
AT1G01720.1,141,870,305,940,950,CDS,GGAGGAGGAGG,G2L1-2
*/
	
	for(int q= 0; q < lines.size(); q++){		

		String line = lines.get(q).trim();

		System.out.println(line);
		String[] tab_split = comma_pattern.split(line);

		if(tab_split.length == 9){
				
			String cur_gene = tab_split[0].trim();
			cur_gene = cur_gene.substring(0,cur_gene.indexOf("."));
			String cur_region = tab_split[6].trim();
			
			if(cur_region.equals("CDS")){
				String cur_type = tab_split[8].trim();
				int type_int = -1;
				
				for(int i = 0; i < gq_type.length; i++){
				
					if(cur_type.equals(gq_type[i])){

						type_int = i;
						break;
					}
				}

				if(type_int == -1 ){
				
				}else{
				
					if(gene_to_k_PDS_dependent_qquad.containsKey(cur_gene)){
				
						Boolean[] exist = gene_to_k_PDS_dependent_qquad.get(cur_gene);
						exist[type_int] = new Boolean(true);
							
						gene_to_k_PDS_dependent_qquad.put(cur_gene,exist);
							
					}else{
						Boolean[] exist = new Boolean[gq_type.length];
						for(int r = 0; r < gq_type.length; r++){
						
							if(r == type_int){
								exist[r] = new Boolean(true);
							}else{

								exist[r] = new Boolean(false);
							}
						}
						gene_to_k_PDS_dependent_qquad.put(cur_gene,exist);
					}
				
				}
					
			}//7
		}//if
			
			
	}//q
	
}//method

//same modules as in medicago_cog_module
static LinkedHashMap<String,String> medi_string_uniprot_to_gff_gene;
static LinkedHashMap<String,String> medi_gff_gene_to_string_uniprot;
static LinkedHashMap<String,String> medi_string_uniprot_to_gff_transcript;

void map_medi_string_uniprot_to_gff_transcript() throws IOException{

	medi_string_uniprot_to_gff_transcript = new LinkedHashMap();
	medi_string_uniprot_to_gff_gene = new LinkedHashMap();
	medi_gff_gene_to_string_uniprot = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("\\|");

	String path = "big/idmapping_2023_09_15.txt";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		if(one_line.contains("ID   ")){

			//System.out.println("498 one_line:" + one_line);
			String uni_id = "";
			String gene_id = "";
			String trans_id = "";
			String gn = "";

			for(int r = q+1; r < lines.size(); r++){
			
				String next_line = lines.get(r).trim();
				
				
				if(next_line.contains("ID   ")){
				
					break;
				}
			
				if(next_line.contains("AC   ")){
				
					uni_id = next_line.replaceFirst("AC   ","").trim();
					uni_id = uni_id.substring(0, uni_id.indexOf(";"));
				}else if(next_line.contains("GN   ")){
				
					
					gn = gn + next_line.substring(next_line.indexOf("   ")+1,next_line.length());
					
					
				}else if(next_line.contains("OS   ")){
				
					String[] str_split = space_pattern.split(gn.trim());
					
					for(int f =0; f < str_split.length; f++){
					
						if(str_split[f].trim().contains("=MTR_")){
						
							gene_id = str_split[f].trim().substring(str_split[f].trim().indexOf("=") +1,str_split[f].trim().length());
							break;
						}
					}
					
					//System.out.println("str_split.length:" + str_split.length);
					for(int f =0; f < str_split.length; f++){
					
						if(str_split[f].trim().contains("EMBL:") || str_split[f].trim().contains("EnsemblPlants:")){
						
							
							String[] bar_split = bar_pattern.split(str_split[f].trim());
							//System.out.println("EMBL:" + bar_split.length);
							
							for(int e = 0; e < bar_split.length; e++){
							
								if(bar_split[e].trim().contains("EMBL:") || bar_split[e].trim().contains("EnsemblPlants:")){
								
									//System.out.println("bar_split[e].trim():" + bar_split[e].trim());
									
									if(bar_split[e].trim().contains(".")){
										trans_id = bar_split[e].trim().substring(bar_split[e].trim().indexOf(":")+1,bar_split[e].trim().indexOf("."));
									}else{
									
										trans_id = bar_split[e].trim().substring(bar_split[e].trim().indexOf(":")+1,bar_split[e].trim().indexOf("}"));
									}
									medi_string_uniprot_to_gff_transcript.put(uni_id, trans_id);
									//System.out.println("error2:" + uni_id + ":" +  trans_id);
									break;
								}
							}
							
							break;
						}
					}
					
					
					
					
					medi_string_uniprot_to_gff_gene.put(uni_id, gene_id);
					//System.out.println("error:" + uni_id + ":" + gene_id);
					medi_gff_gene_to_string_uniprot.put(gene_id,uni_id);
					
					
					break;
				}

			}
		}
	}//lines

				
}//method

static LinkedHashMap<String,String> ara_string_uniprot_to_gff_transcript;

void map_ara_string_uniprot_to_gff_transcript() throws IOException{

	ara_string_uniprot_to_gff_transcript = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("\\|");

	String path = "data/Tair_atxg_nm_uniprot_id_map.idmaps.all";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		String[] space_split = space_pattern.split(one_line);
  		String uni_id = space_split[1].trim();
  		String mn_id = space_split[0].trim();
  		ara_string_uniprot_to_gff_transcript.put(uni_id,mn_id);
			
	}//lines

				
}//method


static LinkedHashMap<String,String> ara_string_uniprot_to_gff_gene;
static LinkedHashMap<String,String> ara_gff_gene_to_string_uniprot;

void map_ara_string_uniprot_to_gff_id() throws IOException{

	ara_string_uniprot_to_gff_gene = new LinkedHashMap();
	ara_gff_gene_to_string_uniprot = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("\\|");

	String path = "data/idmapping_2023_09_26.tsv.id.map";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		String[] space_split = space_pattern.split(one_line);
  		String uni_id = space_split[1].trim();
  		String atxg_id = space_split[0].trim();
		ara_string_uniprot_to_gff_gene.put(uni_id,atxg_id);
		ara_gff_gene_to_string_uniprot.put(atxg_id,uni_id);
		
			
	}//lines

				
}//method

///////////same modules as in medicago_cog_module ends////

static LinkedHashMap<String,Double[]> ara_atxg_to_dna_g4;
static LinkedList<String> ara_g4_x5;

void map_ara_atxg_to_dna_g4() throws IOException{

	ara_atxg_to_dna_g4 = new LinkedHashMap();
	ara_g4_x5 = new LinkedList();
	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/ara.module.e.uniprot.seq.out/score50/";
	File[] files = new File(path).listFiles();
	
	for(int g = 0; g < 2/*files.length*/; g++){
	
		String atxg = files[g].getName();
		atxg = atxg.replaceFirst("ara.module.e.uniprot.all.cog.","");
		atxg = atxg.replaceFirst(".g4","");
		
		List<String> all_lines  = FileUtils.readLines(new File(path + files[g].getName()));
		String one_line = all_lines.get(3).trim();
		String[] split_str = space_pattern.split(one_line);
		int ts = new Integer(split_str[5].trim()).intValue();
		Double[] val = new Double[5];
		val[0] = new Double(ts+0.0);
		
		//1          1582 1587 1600 1605  3  58  GGGTCGGGTGAAAGTGATGGGTCGGG
		
		
		for(int f = 1; f < 4; f++){
		
			val[f] = new Double(split_str[f+1].trim()).doubleValue() - (new Double(split_str[f].trim()).doubleValue()+ts+0.0);
		}
		
		val[4] = new Double(split_str[6].trim());
		
		ara_atxg_to_dna_g4.put(atxg,val);
		System.out.println("atxg g4:"+atxg +":"+ Arrays.toString(val));
		
		if(ts == 5){
		
			ara_g4_x5.add(atxg);
		}
			
		
	}
}//method

static LinkedHashMap<String,Double[]> medi_mtr_to_dna_g4;
static LinkedList<String> medi_g4_x5;

void map_medi_mtr_to_dna_g4() throws IOException{

	medi_mtr_to_dna_g4 = new LinkedHashMap();
	medi_g4_x5 = new LinkedList();
	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/output_medi_module_e/score50/";
	//medi_1064.g4
	
	File[] files = new File(path).listFiles();
	
	for(int g = 0; g < 2/*files.length*/; g++){
	
		String mtr = files[g].getName();
		String input_path = "data/input_medi_module_e/" + mtr.substring(0, mtr.indexOf("."));
		//>MTR_2g045660
		List<String> input_lines  = FileUtils.readLines(new File(input_path));
		String mtr_id = input_lines.get(0);
		mtr_id = mtr_id.replaceFirst(">","");
		
		
		List<String> all_lines  = FileUtils.readLines(new File(path + files[g].getName()));
		String one_line = all_lines.get(3).trim();
		String[] split_str = space_pattern.split(one_line);
		int ts = new Integer(split_str[5].trim()).intValue();
		Double[] val = new Double[5];
		val[0] = new Double(ts+0.0);
		
		//1          1582 1587 1600 1605  3  58  GGGTCGGGTGAAAGTGATGGGTCGGG
		
		for(int f = 1; f < 4; f++){
		
			val[f] = new Double(split_str[f+1].trim()).doubleValue() - (new Double(split_str[f].trim()).doubleValue()+ts+0.0);
		}
		
		val[4] = new Double(split_str[6].trim());
		
		medi_mtr_to_dna_g4.put(mtr_id,val);
		System.out.println("medi g4:"+mtr_id +":"+ Arrays.toString(val));
		
		if(ts == 5){
		
			medi_g4_x5.add(mtr_id);
		}
			
		
	}
}//method



static String[] gquad_range = {"100","200","300","400","500","600","700","800","900","1000","1500","2000","g2000"};
static String[] progs = {"gquad","gquadO","hdna","hdnaO","slipped","str","tfo","zdna"};

static LinkedHashMap<String,Double[][]> ara_atxg_to_gquad_prog;

void map_ara_atxg_to_gquad_prog() throws IOException{

	ara_atxg_to_gquad_prog = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");
	
	String id_path = "data/ara_pfk";
	List<String> ids  = FileUtils.readLines(new File(id_path));
	

	String path = "data/pfk/gquad/ara/";
	List<String>[] all_lines = new List[progs.length];
	
	for(int g = 0; g < progs.length; g++){
	
		all_lines[g]  = FileUtils.readLines(new File(path + "TAIR10_seq_20110103_representative_gene_model_one_line.ipd3.e.module.pfk." + progs[g]));
			
	}
	int[] count = new int[progs.length];
	
	for(int e = 0; e < 2/*ids.size()*/; e++){
	
		String name = ids.get(e);
		Double[][] value = new Double[progs.length][gquad_range.length];
		
		for(int g = 0; g < progs.length; g++){
		
			for(int f = 0; f < gquad_range.length; f++){
				value[g][f] = new Double(0);
			}
		
			for(int f =1 + count[g]; f < all_lines[g].size(); f++){
			
				/*"8"	1	"1886"	"ggtgaagcagagctggattgttttggagaatgcacagtgg"	"40"	"*"
"9"	1	"1981"	"ggatcattcttgttggttaagaggtcaaatcgg"	"33"	"*"
"10"	2	"1204"	"ggcatcggagatgagtcagaaagatccaaactttacatagattggagg"	"48"	"*"
"11"	2	"1447"	"gggtatggttatgggaatttagctggaaggtgttactacacggg"	"44"	"*"

				*/
		
				String one_line = all_lines[g].get(f).trim();
				one_line = one_line.replaceAll("\"","");
				String[] split_str = space_pattern.split(one_line);
				int id_num = new Integer(split_str[1].trim()).intValue();
				
				if(id_num == e+1){
					
					if(!split_str[3].trim().equals("-")){
						
						int width = new Integer(split_str[2].trim()).intValue();
						if(width >= 1 && width < 100){
			
							value[g][0] = value[g][0]+1.0;
						
						}else if(width >= 100 && width < 200){
						
							value[g][1] = value[g][1]+1.0;
						
						}else if(width >= 200 && width < 300){
						
							value[g][2] = value[g][2]+1.0;
						
						}else if(width >= 300 && width < 400){
						
							value[g][3] = value[g][3]+1.0;
						
						}else if(width >= 400 && width < 500){
						
							value[g][4] = value[g][4]+1.0;
						
						}else if(width >= 500 && width < 600){
						
							value[g][5] = value[g][5]+1.0;
						
						}else if(width >= 600 && width < 700){
						
							value[g][6] = value[g][6]+1.0;
						
						}else if(width >= 700 && width < 800){
						
							value[g][7] = value[g][7]+1.0;
						
						}else if(width >= 800 && width < 900){
						
							value[g][8] = value[g][8]+1.0;
						
						}else if(width >= 900 && width < 1000){
						
							value[g][9] = value[g][9]+1.0;
						
						}else if(width >= 1000 && width < 1500){
						
							value[g][10] = value[g][10]+1.0;
						
						}else if(width >= 1500 && width < 2000){
						
							value[g][11] = value[g][11]+1.0;
						
						}else if(width >= 2000){
						
							value[g][12] = value[g][12]+1.0;
						
						}
							
						
					}//if
					
					count[g] = count[g] +1;
				}//if id num
				else if(id_num > e+1){
				
					break;
				}
			}
		}//g
		
		String value_str = "";
		
		for(int i =0; i < value.length; i++){
		
			for(int j = 0; j < value[i].length; j++){
			
				value_str = value_str + value[i][j] + ":" ;
			}
		}
		System.out.println("gquad ara:" + name + ":" + value_str);
		
		ara_atxg_to_gquad_prog.put(name,value);
	}//e
			
	
}//method


static LinkedHashMap<String,Double[][]> medi_mtr_to_gquad_prog;

void map_medi_mtr_to_gquad_prog() throws IOException{

	medi_mtr_to_gquad_prog = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");
	
	String id_path = "data/medi_pfk";
	List<String> ids  = FileUtils.readLines(new File(id_path));
	
	String path = "data/pfk/gquad/medi/";
	List<String>[] all_lines = new List[progs.length];
	
	for(int g = 0; g < progs.length; g++){
	
		all_lines[g]  = FileUtils.readLines(new File(path + "out_get_gene_seq_from_gff_gna.medi.pfk." + progs[g]));
			
	}
	int[] count = new int[progs.length];
	
	for(int e = 0; e < 2/*ids.size()*/; e++){
	
		String name = ids.get(e);
		Double[][] value = new Double[progs.length][gquad_range.length];
		
		for(int g = 0; g < progs.length; g++){
		
			for(int f = 0; f < gquad_range.length; f++){
				value[g][f] = new Double(0);
			}
		
			for(int f =1 + count[g]; f < all_lines[g].size(); f++){
			
				/*"8"	1	"1886"	"ggtgaagcagagctggattgttttggagaatgcacagtgg"	"40"	"*"
"9"	1	"1981"	"ggatcattcttgttggttaagaggtcaaatcgg"	"33"	"*"
"10"	2	"1204"	"ggcatcggagatgagtcagaaagatccaaactttacatagattggagg"	"48"	"*"
"11"	2	"1447"	"gggtatggttatgggaatttagctggaaggtgttactacacggg"	"44"	"*"

				*/
		
				String one_line = all_lines[g].get(f).trim();
				one_line = one_line.replaceAll("\"","");
				String[] split_str = space_pattern.split(one_line);
				int id_num = new Integer(split_str[1].trim()).intValue();
				
				if(id_num == e+1){
					
					if(!split_str[3].trim().equals("-")){
						
						int width = new Integer(split_str[2].trim()).intValue();
						if(width >= 1 && width < 100){
			
							value[g][0] = value[g][0]+1.0;
						
						}else if(width >= 100 && width < 200){
						
							value[g][1] = value[g][1]+1.0;
						
						}else if(width >= 200 && width < 300){
						
							value[g][2] = value[g][2]+1.0;
						
						}else if(width >= 300 && width < 400){
						
							value[g][3] = value[g][3]+1.0;
						
						}else if(width >= 400 && width < 500){
						
							value[g][4] = value[g][4]+1.0;
						
						}else if(width >= 500 && width < 600){
						
							value[g][5] = value[g][5]+1.0;
						
						}else if(width >= 600 && width < 700){
						
							value[g][6] = value[g][6]+1.0;
						
						}else if(width >= 700 && width < 800){
						
							value[g][7] = value[g][7]+1.0;
						
						}else if(width >= 800 && width < 900){
						
							value[g][8] = value[g][8]+1.0;
						
						}else if(width >= 900 && width < 1000){
						
							value[g][9] = value[g][9]+1.0;
						
						}else if(width >= 1000 && width < 1500){
						
							value[g][10] = value[g][10]+1.0;
						
						}else if(width >= 1500 && width < 2000){
						
							value[g][11] = value[g][11]+1.0;
						
						}else if(width >= 2000){
						
							value[g][12] = value[g][12]+1.0;
						
						}
							
						
					}//if
					
					count[g] = count[g] +1;
				}//if id num
				else if(id_num > e+1){
				
					break;
				}
			}
		}//g
		
		String value_str = "";
		
		for(int i =0; i < value.length; i++){
		
			for(int j = 0; j < value[i].length; j++){
			
				value_str = value_str + value[i][j] + ":" ;
			}
		}
		//System.out.println("gquad:" + name + ":" + value_str);
		
		System.out.println("medi gquad:" + name + ":" + value_str);
		
		medi_mtr_to_gquad_prog.put(name,value);
	}//e
			
	
}//method

static String[] prot_metal = {"MetalBinding","Ca","Co","Cu","Fe","K","Mg","Mn","Na","Ni","Zn"};


static LinkedHashMap<String,Double[]> ara_atxg_to_prot_metal_binding;

void map_ara_atxg_to_prot_metal_binding() throws IOException{

	ara_atxg_to_prot_metal_binding = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/mymetal/output_ara/";
	//ara.module.e.uniprot.all.cog.AT1G01180.out
	
	File[] files = new File(path).listFiles();
	
	for(int g = 0; g < 2/*files.length*/; g++){
	
		String atxg_id = files[g].getName();
		atxg_id = atxg_id.replaceFirst("ara.module.e.uniprot.all.cog.","");
		atxg_id = atxg_id.replaceFirst(".out","");
		
		
		List<String> all_lines  = FileUtils.readLines(new File(path + files[g].getName()));
		String one_line = all_lines.get(8).trim();
		String[] split_str = space_pattern.split(one_line);
		Double[] val = new Double[prot_metal.length];
			
		//3880.G7IT21	0.13	0.14	0.08	0.0	0.0	0.04	0.23	0.01	0.03	0.22	0.06
			
			
		for(int f = 1; f < split_str.length; f++){
			
			val[f-1] = new Double(split_str[f].trim());
		}
			
		ara_atxg_to_prot_metal_binding.put(atxg_id,val);
		System.out.println("ara prot metal:"+atxg_id +":"+ Arrays.toString(val));
		
	}
}//method

static LinkedHashMap<String,Double[]> medi_mtr_to_prot_metal_binding;

void map_medi_mtr_to_prot_metal_binding() throws IOException{

	medi_mtr_to_prot_metal_binding = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/mymetal/output_medi/";
	//medi.module.e.uniprot.all.cog.S5MFQ5.out
	
	File[] files = new File(path).listFiles();
	
	for(int g = 0; g < 2/*files.length*/; g++){
	
		String uni_id = files[g].getName();
		uni_id = uni_id.replaceFirst("medi.module.e.uniprot.all.cog.","");
		uni_id = uni_id.replaceFirst(".out","");
		
		String mtr_id = medi_string_uniprot_to_gff_gene.get(uni_id);
		
		if(mtr_id != null){
			List<String> all_lines  = FileUtils.readLines(new File(path + files[g].getName()));
			String one_line = all_lines.get(8).trim();
			String[] split_str = space_pattern.split(one_line);
			Double[] val = new Double[prot_metal.length];
			
			//3880.G7IT21	0.13	0.14	0.08	0.0	0.0	0.04	0.23	0.01	0.03	0.22	0.06
			
			
			for(int f = 1; f < split_str.length; f++){
			
				val[f-1] = new Double(split_str[f].trim());
			}
			
			
			medi_mtr_to_prot_metal_binding.put(mtr_id,val);
			System.out.println("medi prot metal:"+mtr_id +":"+ Arrays.toString(val));
			
		}
			
		
	}
}//method

}//inner class. dna structure and prot metal binding

static class ddi_interface_anno_module{


//some ID mapping functions (e.g., pdbid to pfam, atxg id to pfam etc.) were from art_founddation_log program

static String[] hmm_aa = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};

static String[] hmm_indel = {"m->m", "m->i" , "m->d" ,"i->m" ,  "i->i"  , "d->m" ,  "d->d"};
	
static String[] gcg_header = {"Cons", "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","U","B","Z","X",   "Gap","Len"};

static LinkedHashMap<String,String> domain_to_pfam;
static LinkedHashMap<String,String> pfam_to_domain;
void map_domain_pfam_3did() throws IOException{

	Pattern space_pattern = Pattern.compile("\\s+");

	domain_to_pfam = new LinkedHashMap();
	pfam_to_domain = new LinkedHashMap();

	String path = "data/map_domain_pfam.final.2";

	List<String> lines  = FileUtils.readLines(new File(path));

	for(int q= 0; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();
		
		String[] split_str = space_pattern.split(one_line);
		String name = split_str[0].trim();
		String pfam = split_str[1].trim();
		pfam = pfam.substring(0,pfam.indexOf("."));
		
		domain_to_pfam.put(name,pfam);
		pfam_to_domain.put(pfam,name);
	}//q

	//test

	/*List<String> keys_list = new LinkedList<String>(domain_to_pfam.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		System.out.println("domain_to_pfam:" + keys_list.get(i) + ":" + domain_to_pfam.get(keys_list.get(i)));
	}

	keys_list = new LinkedList<String>(pfam_to_domain.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		System.out.println("pfam_to_domain:" + keys_list.get(i) + ":" + pfam_to_domain.get(keys_list.get(i)));
	}*/
}

//argument:domain dimer ex.120_Rick_ant@120_Rick_ant
//return: how many topology this dimer contains

static int map_interface_topo_num(String filename) throws IOException{

	if(new File("data/3did/interface_flat/" + filename + ".txt").exists()){
	    	List<String> lines  = FileUtils.readLines(new File("data/3did/interface_flat/" + filename + ".txt"));
	    
		//System.out.println("topo num error:"+"/media/ubuntu/frontiers256/frontiers/3did/interface_flat/" + filename + ".txt" + ":" + (lines.size()-2));
		String last = lines.get(lines.size()-2).trim();
		String num = last.substring(0,last.indexOf(":"));
		num = num.replaceFirst("#=IF","").trim();
		
		int return_val = new Integer(num).intValue();
		return(return_val);
	}
	
	return(0);

}

static String[] get_cons_hmm_GCG_profile(String pfam_id) throws IOException{

	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");

	//input directory is too big. can be made with .hmm file downloaded from pfam website. covert hmm to prf using cran R package. only one pfam file is given as an example
	//String path = "data/pfam_hmm/indi_file/" + pfam_id + ".txt.prf";
	String path = "/media/ubuntu/frontiers256/frontiers/pfam/pfam_hmm/indi_file/" + pfam_id + ".txt.prf";
	/*
!!AA_PROFILE 1.0
(Peptide) HMMCONVERT v2.3.2 Length: 264 7tm_1|PF00001.24|7 transmembrane receptor (rhodopsin family)
   Profile converted from a profile HMM using HMMER v2.3.2 emulation.
   Use -nonor -noave -gap=10 -len=1 with profilesearch and friends
      to get the closest approximation to HMMER bit scores.
   WARNING: There is a loss of information in this conversion.
      Neither the scores nor even the rank order of hits will be precisely
      preserved in a comparison of HMMER hmmsearch to GCG profilesearch.
      The profile score is an approximation of the (single-hit) HMMER score.

Cons    A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y     U     B     Z     X   Gap   Len ..
 G    -19  -134   -61    43  -171   223   -36  -144    -2  -140   -74     5  -166     6   -54     4   -26  -101   -42  -124     4   -32    29   -42   100   100
 N    -83  -191   -95  -131  -332  -186  -193  -345  -193  -363  -292   415  -245  -167  -232  -130  -154  -275  -337  -285  -130   127  -144  -189    70   111
 L     -5    48  -276  -220    -3   -30   -98   108  -175   137    57  -155  -211  -138  -167   -37    14   106   -73   -44   -37  -223  -190   -51    71   111
 
*/
    
  
		
    List<String> lines  = FileUtils.readLines(new File(path));
    LinkedList<String> input = new LinkedList();
	//G    -19  -134   -61    43  -171   223   -36  -144    -2  -140   -74     5  -166     6   -54     4   -26  -101   -42  -124     4   -32    29   -42   100   100

	
	for(int q= 0; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();

		if(one_line.startsWith("Cons")){

			for(int w= q+1; w < lines.size(); w++){
	
				String one_line2 = lines.get(w).trim();
				String[] space_split = space_pattern.split(one_line2);
				String cons = space_split[0].trim();
				input.add(cons);
			}
			
		}
		
			
 	}//q

	String[] ret_val = input.toArray(new String[0]);

	

	return(ret_val);
 

}


static String[] aa_list = {"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"};
//among hetero dom dimer (e.g., dom1dom2), the max topology number is 117. in other words, there are 117 IF lines for an interface flat file.
//argument:domain dimer ex.120_Rick_ant@120_Rick_ant
//return: per aa, sum residue of hmm of the first domain,how many topology  in all domain dimer interfaces. 
//if this residue belongs to a topology, frequency goes up.

static int[] map_inter_topo(String filename) throws IOException{

	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern dash_pattern = Pattern.compile("-");

	System.out.println("map_inter_topo:"+filename);
	int[] return_val = new int[1];
	if(new File("data/3did/interface_flat/" + filename + ".txt").exists()){

		List<String> lines  = FileUtils.readLines(new File("data/3did/interface_flat/" + filename + ".txt"));

		System.out.println("map_inter_topo:"+"data/3did/interface_flat/" + filename + ".txt");
		
		//get how many seq in hmm profile
		String dom1 = filename.substring(0,filename.indexOf("@"));
		String dom2 = filename.substring(filename.indexOf("@")+1,filename.length());
		
		if(!dom1.equals(dom2)){
			String pfam1 = domain_to_pfam.get(dom1);
			
			System.out.println("pfam1:" + pfam1 + ":" + dom1);
			String[] cons = get_cons_hmm_GCG_profile(pfam1);

			System.out.println("cons size:" + cons.length);

			//get topo_num
			int topo_num = map_interface_topo_num(filename);
			System.out.println("topo_num:" + topo_num);

			return_val = new int[aa_list.length];
			boolean[][] temp = new boolean[aa_list.length][topo_num+1];
			
			/*
			#=ID	ZZ	Myb_DNA-binding
		#=IF	0: 1-1 3-3 5-5 19-19 21-21 23-23 26-27
		#=IF	1: 5-5 15-19 21-21 23-23 26-27

		#=ID	DNApol3-delta_C	DNApol3-delta_C
		#=IF	0: 54-54 74-81
		//
		//

			*/

			int max = 0;
			for(int q= 1; q < lines.size()-1; q++){
			
				String one_line = lines.get(q).trim();
				
				
				if(one_line.contains("IF")){

					one_line = one_line.replaceFirst("#=IF","").trim();
					String[] split_colon = colon_pattern.split(one_line);
					int topo_indi = new Integer(split_colon[0].trim()).intValue() -1;
					System.out.println("topo_indi:" + topo_indi);
					
					String[] split_str = space_pattern.split(split_colon[1].trim());

					for(int w = 0; w < split_str.length; w++){

						
						String[] split_range = dash_pattern.split(split_str[w].trim());
						int start = new Integer(split_range[0].trim()).intValue();
						int end = new Integer(split_range[1].trim()).intValue();

						if(max < end){
							max = end;
						}
						for(int e = start; e < end+1; e++){

							int cons_num = e-1;
							String cons_indi = cons[cons_num];

							int aa_int = -1;
							for(int r = 0; r < aa_list.length; r++){

								if(cons_indi.equals(aa_list[r])){

									aa_int = r;
									break;
								}
							}
							
							//System.out.println("error:" + "data/3did/interface_flat/" + filename + ".txt" + ":" + cons_indi  + ":" +cons_num+":"+ aa_int + ":" + topo_indi);
							if(aa_int != -1 && topo_indi != -1){
								temp[aa_int][topo_indi] =true;
								System.out.println("return_val:" + e + ":" + topo_indi + ":true" );
							}

						}

					}

					
				}//if
				
		 	}//q

			System.out.println("max:" + filename + ":" + max);
			for(int q = 0; q < temp.length; q++){

				int count = 0;
				for(int w = 0; w < temp[q].length; w++){

					if(temp[q][w]){

						count = count +1;
					}
				}
				
				return_val[q] = count;
			}
		}//different
	}//exist
	 
	return(return_val);
}

//string filename is from protcad nr pfam.
//need to decide which pfam should be the one that needs hmm information
//then, for example hmm from TF and other BP from LLPS factor and so on..
//let us assume that all TF has LLPS properties
//if TF LLPS property is weak, it may associate with strong LLPS factor, and so on

//if TF has no LLPS property, then what are the different behaviors?
//in terms of histone protein modification, are there any differences in target behaviors?
//LLPS factor associate with some proteins that interact with epigenetic factors that affect histone proteins? -- these need to be answered by literature search. no time.. just make modules and try later with biomarkers

static LinkedHashMap<String, Boolean[]> only_onebp;

static LinkedHashMap<LinkedList<String>, Boolean[]> only_multibp;

static LinkedHashMap<String, LinkedList<LinkedList<String>>> pfam_to_sub_folder;

static LinkedHashMap<LinkedList<String>, Boolean[]> one_multibp;

static String[] first_cat_list = {"onebp","multiplebp","bothbp"};
static int[] first_cat_num_list = {0,1,2};

static String[] cat_stat = {"non_nr_dmi_nr","nr_dmi_nr","non_nr_dmi_non_nr","nr_dmi_non_nr","non_nr","nr"};
//first refers properties of ddi, for example, non_nr_dmi_nr. first part (non_nr) refers to ddi and the second part- dmi_nr refers that dmi is nr. if no dmi, then, the last two values of the array ("non_nr", "nr" are used.

void map_cat_stat() throws IOException{

	//nr_domain: dom1@dom2@dom3
	Pattern at_pattern = Pattern.compile("@");
	only_onebp = new LinkedHashMap();
	only_multibp =  new LinkedHashMap();
	one_multibp =  new LinkedHashMap();
	pfam_to_sub_folder  =  new LinkedHashMap();
	
	String dir = "data/3did/out_global_organize";
	List<String> dir_lines  = FileUtils.readLines(new File(dir));

	for(int z = 0; z < dir_lines.size(); z++){

		String pfam_name = dir_lines.get(z);
		
		if(new File("big/global_inter_organize/" + pfam_name).exists()){
		
			//System.out.println("map:exist:"+pfam_name);
			//System.out.println("map:exist domain:"+pfam_name +":"+ pfam_to_domain.get(pfam_name));
		
			File[] multi_onebp = new File("big/global_inter_organize/" + pfam_name).listFiles();

			//System.out.println("map:exist:multi_onebp:"+multi_onebp.length);
			
			int first_cat = -1;
			
			LinkedList<LinkedList<String>> all_sub_folder= new LinkedList();
			
			if(multi_onebp.length == 1){

				int start = 0;

				if(multi_onebp[0].getName().equals("onebp")){

					first_cat = 0;
					

				}else{
					first_cat = 1;
					

				}

				File[] bp_files = new File("big/global_inter_organize/" + pfam_name + "/" + first_cat_list[first_cat] + "/ddi/").listFiles();
				//System.out.println("map:exist:bp_files:"+"/media/ubuntu/frontiers256/frontiers/3did/global_inter_organize/" + pfam_name + "/" + first_cat_list[first_cat] + "/ddi/");
			
			
				for(int a = 0; a < bp_files.length; a++){
				
					LinkedList<String> name_list = new LinkedList();
					name_list.add(pfam_name);
					
					if(bp_files[a].getName().contains("@")){
						String[] key_split = at_pattern.split(bp_files[a].getName());
						for(int b = 0; b < key_split.length; b++){

							name_list.add(key_split[b].trim());
						}
						
					}
					
					all_sub_folder.add(name_list);
					File[] nr_folder = bp_files[a].listFiles();
					
					//System.out.println("map:exist:nr_folder:"+nr_folder.length);
					
					Boolean[] cat_val = new Boolean[cat_stat.length];
					for(int c = 0; c < cat_stat.length; c++){

						cat_val[c]= new Boolean(false);
					}

					for(int b = 0; b < nr_folder.length; b++){

						for(int c = 0; c < cat_stat.length; c++){
			

							if(nr_folder[b].getName().equals(cat_stat[c])){
								cat_val[c]= new Boolean(true);
								break;
							}
						}
					}
					//System.out.println("map:exist:cat_val:"+Arrays.toString(cat_val));
					if(first_cat == 0){
						only_onebp.put(pfam_name,cat_val);
						
						if(only_onebp.containsKey(pfam_name)){

							Boolean[] exist = only_onebp.get(pfam_name);
							
							for(int p =0; p < exist.length; p++){
							
								if(cat_val[p]){
									exist[p] = true;
								}
							}											
							only_onebp.put(pfam_name,exist);
						}else{
							only_onebp.put(pfam_name,cat_val);
						}

					}else{

						Collections.sort((List)name_list);
						if(only_multibp.containsKey(name_list)){

							Boolean[] exist = only_multibp.get(name_list);
							
							for(int p =0; p < exist.length; p++){
							
								if(cat_val[p]){
									exist[p] = true;
								}
							}											
							only_multibp.put(name_list,exist);
						}else{
							only_multibp.put(name_list,cat_val);
						}
					}
				}

				//
				pfam_to_sub_folder.put(pfam_name,all_sub_folder);

			}else if(multi_onebp.length == 2){

				first_cat = 2;
				

				//for(int p = 0; p < 2; p++){

					File[] bp_files = new File("big/global_inter_organize/" + pfam_name + "/" + first_cat_list[1] + "/ddi/").listFiles();
					//System.out.println("/media/ubuntu/frontiers256/frontiers/3did/global_inter_organize/" + pfam_name + "/" + first_cat_list[1] + "/ddi/");
					for(int a = 0; a < bp_files.length; a++){
					
						LinkedList<String> name_list = new LinkedList();
						name_list.add(pfam_name);
						
						if(bp_files[a].getName().contains("@")){
							String[] key_split = at_pattern.split(bp_files[a].getName());
							for(int b = 0; b < key_split.length; b++){

								name_list.add(key_split[b].trim());
							}
							
						}

						all_sub_folder.add(name_list);
						
						File[] nr_folder = bp_files[a].listFiles();
						Boolean[] cat_val = new Boolean[cat_stat.length];
						
						//System.out.println("multi map:exist:nr_folder:"+nr_folder.length);
						
						for(int c = 0; c < cat_stat.length; c++){

							cat_val[c]= new Boolean(false);
						}

						
						for(int b = 0; b < nr_folder.length; b++){

							for(int c = 0; c < cat_stat.length; c++){
				

								if(nr_folder[b].getName().equals(cat_stat[c])){
									cat_val[c]= new Boolean(true);
									break;
								}
							}
						}//b

						//System.out.println("multi map:exist:cat_val:"+Arrays.toString(cat_val));
						
						Collections.sort((List)name_list);
						if(one_multibp.containsKey(name_list)){

							Boolean[] exist = one_multibp.get(name_list);
							
							for(int p =0; p < exist.length; p++){
							
								if(cat_val[p]){
									exist[p] = true;
								}
							}											
							one_multibp.put(name_list,exist);
						}else{
							one_multibp.put(name_list,cat_val);
						}
					}
					
					//
				pfam_to_sub_folder.put(pfam_name,all_sub_folder);
				//}//p
			}//if
		}
	}//z

	//test
	/*
	List<String> keys_list = new LinkedList<String>(only_onebp.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		
		//System.out.println("only_onebp:" + keys_list.get(i) + ":" + Arrays.toString(ArrayUtils.toPrimitive(only_onebp.get(keys_list.get(i)))));
	}

	keys_list = new LinkedList<String>(only_multibp.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		
		//System.out.println("only_multibp:" + keys_list.get(i) + ":" + Arrays.toString(ArrayUtils.toPrimitive(only_multibp.get(keys_list.get(i)))));
	}

	keys_list = new LinkedList<String>(one_multibp.keySet());
	////System.out.println("one_multibp:error:" + keys_list.size());

	for(int i = 0; i < keys_list.size(); i++){

		//System.out.println("one_multibp:" + keys_list.get(i) + ":" + Arrays.toString(ArrayUtils.toPrimitive(one_multibp.get(keys_list.get(i)))));
	}
	*/
	
}//method


//domains: pfam_name
static boolean[] get_cat_stat_only_onebp(String pfam_name) throws IOException{

	//Pattern at_pattern = Pattern.compile("@");

	//LinkedList<String> query = new LinkedList();
	//String[] key_split = at_pattern.split(domains);
	//for(int q = 0; q < key_split.length; q++){

		//query.add(key_split[q].trim());
	//}
	Boolean[] return_obj = null;

	List<String> keys_list = new LinkedList<String>(only_onebp.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		String key = keys_list.get(i);
		//key_split = at_pattern.split(key);
		//LinkedList<String> key_name = new LinkedList();
		//for(int q = 0; q < key_split.length; q++){

			//key_name.add(key_split[q].trim());
		//}
		
		//if(key_name.containsAll(query) && query.containsAll(key_name) && query.size() == key_name.size()){
		
		if(key.equals(pfam_name)){

			return_obj = only_onebp.get(key);
			break;
		}
	}//i

	boolean[] return_val = null;
	if(return_obj != null){
		return_val = ArrayUtils.toPrimitive(return_obj);
	}
	return(return_val);

}

//domains: pfam_name + "@" + bp_files[a].getName()
static boolean[] get_cat_stat_only_multibp(String domains) throws IOException{

	Pattern at_pattern = Pattern.compile("@");

	LinkedList<String> query = new LinkedList();
	String[] key_split = at_pattern.split(domains);
	for(int q = 0; q < key_split.length; q++){

		query.add(key_split[q].trim());
	}
	Collections.sort((List)query);
	Boolean[] return_obj = only_multibp.get(query);

	/*List<String> keys_list = new LinkedList<String>(only_multibp.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		String key = keys_list.get(i);
		key_split = at_pattern.split(key);
		LinkedList<String> key_name = new LinkedList();
		for(int q = 0; q < key_split.length; q++){

			key_name.add(key_split[q].trim());
		}
		
		if(key_name.containsAll(query) && query.containsAll(key_name) && query.size() == key_name.size()){

			return_obj = only_multibp.get(key);
			break;
		}
	}//i
	*/

	boolean[] return_val = null;
	if(return_obj != null){
		return_val = ArrayUtils.toPrimitive(return_obj);
	}

	return(return_val);

}

//domains: pfam_name + "@" + bp_files[a].getName()
static boolean[] get_cat_stat_one_multibp(String domains) throws IOException{

	Pattern at_pattern = Pattern.compile("@");

	LinkedList<String> query = new LinkedList();
	String[] key_split = at_pattern.split(domains);
	for(int q = 0; q < key_split.length; q++){

		query.add(key_split[q].trim());
	}
	Collections.sort((List)query);
	Boolean[] return_obj = one_multibp.get(query);

	/*List<String> keys_list = new LinkedList<String>(one_multibp.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		String key = keys_list.get(i);
		key_split = at_pattern.split(key);
		LinkedList<String> key_name = new LinkedList();
		for(int q = 0; q < key_split.length; q++){

			key_name.add(key_split[q].trim());
		}
		
		if(key_name.containsAll(query) && query.containsAll(key_name) && query.size() == key_name.size()){

			return_obj = one_multibp.get(key);
			break;
		}
	}//i
	*/

	boolean[] return_val = null;
	if(return_obj != null){
		return_val = ArrayUtils.toPrimitive(return_obj);
	}

	return(return_val);

}

//argument: domains (main pfam + domains) same as get_cat_stat_one_multibp
//return: whether this nr domain has onebp, multipb, or both
static int getFirstCatNum (String domains) throws IOException{

	
	boolean[] cat_stat_cur = get_cat_stat_only_onebp(domains);
	int cat_stat_int= -1;
	
	if(cat_stat_cur == null){

		cat_stat_cur = get_cat_stat_only_multibp(domains);

		if(cat_stat_cur == null){

			cat_stat_cur = get_cat_stat_one_multibp(domains);
			if(cat_stat_cur != null){
				cat_stat_int= 2;

			}

		}else{

			cat_stat_int= 1;

		}

	}else{
		cat_stat_int= 0;
	}

	return(cat_stat_int);
}

//argument: domains (main pfam + domains) same as get_cat_stat_one_multibp
//return: whether this domain folder has the following subfolders- "non_nr_dmi_nr","nr_dmi_nr","non_nr_dmi_non_nr","nr_dmi_non_nr","non_nr","nr"
//if it has the folder, returns true
static boolean[] getNrFolder (int cat_stat_int, String domains) throws IOException{

	boolean[] cat_stat_cur = null;
	if(cat_stat_int == 0){
		cat_stat_cur = get_cat_stat_only_onebp(domains);

	}else if(cat_stat_int == 1){

		cat_stat_cur = get_cat_stat_only_multibp(domains);
	}else if(cat_stat_int == 2){
		cat_stat_cur = get_cat_stat_one_multibp(domains);
	}
	return(cat_stat_cur);
}

		
//among hetero dom dimer (e.g., dom1dom2), the max topology number is 117. in other words, there are 117 IF lines for an interface flat file.
//argument:domain dimer ex.120_Rick_ant@120_Rick_ant
//return: per aa, sum residue of hmm of the first domain,how many topology  in all domain dimer interfaces. 
//if this residue belongs to a topology, frequency goes up.

//argument: domains (nr_pfam + domains) same as get_cat_stat_one_multibp
//return: double[][][] (row: residue num, col: topology num*domain number) per folder ("non_nr_dmi_nr","nr_dmi_nr","non_nr_dmi_non_nr","nr_dmi_non_nr","non_nr","nr"): first dim: folder num, second dim: aa. if(for topology that has hit in BD section, aa has average  (PS value), otherwise, 0.0). if one_multibp, only multibp returns
		
static double[][] getOrderedDomDimerGlobalInterfacePS(String domains) throws IOException{

	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern dash_pattern = Pattern.compile("-");
	Pattern at_pattern = Pattern.compile("@");
	Pattern comma_pattern = Pattern.compile(",");

	String pfam1 = domains.substring(0,domains.indexOf("@"));
	//System.out.println("pfam1:" + pfam1);
	
	String bps = domains.substring(domains.indexOf("@")+1,domains.length());
	//System.out.println("bps:" + bps);

	String[] cons =get_cons_hmm_GCG_profile(pfam1);
	
	int residue_num = cons.length;
	//System.out.println("residue_num:" + residue_num);
	
	int cat_stat_int= getFirstCatNum(domains);
	//System.out.println("cat_stat_int:" + cat_stat_int);
	boolean[] cat_stat_cur = getNrFolder(cat_stat_int,domains);
	//System.out.println("cat_stat_cur:" + cat_stat_cur[0] + ":" + cat_stat_cur[1] + ":" +cat_stat_cur[2] + ":" +cat_stat_cur[3] + ":" +cat_stat_cur[4] + ":" +cat_stat_cur[5]);

	LinkedList<Double>[][] temp_val = new LinkedList[6][];
	double[][] return_val = new double[1][1];
	if((cat_stat_int ==1 || cat_stat_int ==2) && cat_stat_cur != null){

		String dir = "big/global_inter_organize/";
		////get global min per bps folder. later calcul average + min
		String main_dom = pfam_to_domain.get(pfam1);
		
		String path = dir + pfam1 + "/" + first_cat_list[1] + "/ddi/" + bps;
		//System.out.println("path:" +path);
		
		if(new File(path).exists()){
			
			for(int i =0;i < cat_stat.length; i++){


				if(cat_stat_cur[i]){

					if(first_cat_list[cat_stat_int].equals("bothbp")){

						first_cat_list[cat_stat_int] = "multiplebp";
					}

					path = dir + pfam1 + "/" + first_cat_list[cat_stat_int] + "/ddi/" + bps + "/" + cat_stat[i];

					//System.out.println("if(cat_stat_cur[i]) path:" +path);

					File[] bp_sub_files = new File(path).listFiles();
					
					temp_val[i] = new LinkedList[aa_list.length];
					for(int a = 0; a < aa_list.length; a++){
						temp_val[i][a] = new LinkedList();

					}

					

					for(int a = 0; a < bp_sub_files.length; a++){

						//readlines
						List<String> lines  = FileUtils.readLines(bp_sub_files[a]);
						

							String ps_line = lines.get(2);
							ps_line = ps_line.replaceFirst("#=PS","").trim();
							//System.out.println("ps_line:" +ps_line);
							String[] split_ps = comma_pattern.split(ps_line);
							double[] ps = new double[residue_num];
							for(int w = 0; w < split_ps.length; w++){

								String indi = split_ps[w].trim();
								String[] split_indi = space_pattern.split(indi);
								int cons_num = new Integer(split_indi[0].trim()).intValue()-1;

								
								String cons_indi = cons[cons_num];
								//System.out.println(cons_indi);

								int aa_int = -1;
								for(int r = 0; r < aa_list.length; r++){

									if(cons_indi.equals(aa_list[r])){

										aa_int = r;
										break;
									}
								}

								if(aa_int != -1){
									String val = split_indi[1].trim();
									val = val.replaceFirst("\\(","");
									val = val.replaceFirst("\\)","").trim();
									Double hmm_ps = new Double(val);
									temp_val[i][aa_int].add(hmm_ps);
									//System.out.println("hmm_ps:" +hmm_ps + ":" + cons_num);
									//11 (0.2)
								}

							}

						//per residue per topo per domain dimer --> add values to linkedlist<double>, later get average + min
					}//a

				}// if(cat_stat_cur[i]){

				
			}//i	

			//now.. get average. prepare return values
			
			
			
			return_val = new double[6][aa_list.length];	
			for(int i =0;i < temp_val.length; i++){

				if(temp_val[i] != null){

					for(int g = 0; g < temp_val[i].length; g++){

						if(temp_val[i][g].size() == 0){

							return_val[i][g] = 0.0;
						}else{
						
							double[] values =  ArrayUtils.toPrimitive(temp_val[i][g].toArray(new Double[0]));
						
							if(values.length == 1){

								return_val[i][g] = values[0];
								 												
								
							}else if(values.length >1){
						
		
								DescriptiveStatistics desc = new DescriptiveStatistics(values);

								return_val[i][g] = desc.getMean();
								//System.out.println(">1:" +return_val[i][g]);
								
							}
							
						}
					}//g
				}//not null
			}//i

			return(return_val);

		}//exist
	}//if 1 or 2
	
	return(return_val);
}//method	



static LinkedList<Double[]> map_hmm_profile(String pfam_id) throws IOException{

	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");


	String path = "data/pfam_hmm/indi_file/" + pfam_id + ".txt.hmm";
	/*
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   4.99131  6.56085  5.46586  2.59346  4.30426  5.52625  3.96988  3.30179  1.66209  3.89165  5.93737  4.48998  0.99769  4.28392  5.00216  3.69251  4.47318  2.42411  6.93622  2.29751
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.10230  6.13491  2.35305  0.61958  0.77255  0.00000        *
      1   6.63533  7.64427  7.13096  7.21415  7.80411  6.25106  7.83759  8.12206  7.42137  7.28282  8.55037  7.42568  0.01255  7.71438  7.24074  6.92562  7.12904  7.60682  8.22979  8.01779      1 P - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00356  6.03616  6.75851  0.61958  0.77255  0.23681  1.55658
......
     10   5.05052  7.12404  5.02035  4.50912  6.82210  5.38727  5.23333  6.05134  0.33891  5.31409  6.23416  2.50813  5.77649  2.21185  3.52082  4.97419  5.14638  5.73888  7.21478  6.15950     10 K - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.00217  6.13386        *  0.61958  0.77255  0.00000        *
*/
    
  
		
    List<String> lines  = FileUtils.readLines(new File(path));
    LinkedList<Double[]> input = new LinkedList();
    
	for(int q= 27; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();
		//System.out.println("one_line:" + one_line);
		String next_line = "";
		String next_line2 = "";
		
		/*if(q == 24){

			one_line = one_line.replaceFirst("COMPO","").trim();
		}*/

		
		String[] split_str = space_pattern.split(one_line);

		if(split_str.length > 20){

			//System.out.println("split_str.length > 20:" + one_line);
			Double[] arr = new Double[47];

			for(int w = 1; w < 21; w++){

				arr[w-1] = new Double(split_str[w].trim());

			}

			if(q < lines.size()-3){
				next_line = lines.get(q+1).trim();
				next_line2 = lines.get(q+2).trim();
				//System.out.println("q < lines.size()-3 next_line:" + next_line);
				//System.out.println("q < lines.size()-3 next_line2:" + next_line2);

			}

			String[] split_next = space_pattern.split(next_line);

			if(split_next.length == 20){

				//System.out.println("split_next.length == 20 next_line:" + next_line);
				for(int w = 0; w < split_next.length; w++){

					arr[w+20] = new Double(split_next[w].trim());

				}
			}

			String[] split_next2 = space_pattern.split(next_line2);
			if(split_next2.length == 7){
				
				//System.out.println("split_next2.length == 7 next_line2:" + next_line2);

				for(int w = 0; w < split_next2.length; w++){

					if(split_next2[w].trim().equals("*")){

						arr[w+40] = Double.NaN;
					}else{

						arr[w+40] = new Double(split_next2[w].trim());
					}

				}
			}
			input.add(arr);
		}//20
		
		
			
 	}//q

	/*for(int j = 0; j < input.size(); j++){

		Double[] arr = input.get(j);

		for(int i = 0; i < arr.length; i++){
			//System.out.println("final:" + j + ":" + i + ":" + arr[i]);
		}
	}*/

	return(input);
 

}

static LinkedList<Double[]> map_hmm_GCG_profile(String pfam_id) throws IOException{

	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");


	String path = "data/pfam_hmm/indi_file/" + pfam_id + ".txt.prf";
	/*
!!AA_PROFILE 1.0
(Peptide) HMMCONVERT v2.3.2 Length: 264 7tm_1|PF00001.24|7 transmembrane receptor (rhodopsin family)
   Profile converted from a profile HMM using HMMER v2.3.2 emulation.
   Use -nonor -noave -gap=10 -len=1 with profilesearch and friends
      to get the closest approximation to HMMER bit scores.
   WARNING: There is a loss of information in this conversion.
      Neither the scores nor even the rank order of hits will be precisely
      preserved in a comparison of HMMER hmmsearch to GCG profilesearch.
      The profile score is an approximation of the (single-hit) HMMER score.

Cons    A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y     U     B     Z     X   Gap   Len ..
 G    -19  -134   -61    43  -171   223   -36  -144    -2  -140   -74     5  -166     6   -54     4   -26  -101   -42  -124     4   -32    29   -42   100   100
 N    -83  -191   -95  -131  -332  -186  -193  -345  -193  -363  -292   415  -245  -167  -232  -130  -154  -275  -337  -285  -130   127  -144  -189    70   111
 L     -5    48  -276  -220    -3   -30   -98   108  -175   137    57  -155  -211  -138  -167   -37    14   106   -73   -44   -37  -223  -190   -51    71   111
...
! 261
 P   -347  -310  -390  -428  -458  -332  -392  -521  -437  -501  -468  -384   427  -425  -421  -359  -367  -465  -380  -444  -359  -387  -427  -376    71   111
 I    -75   -34  -302  -246   122  -187  -115   235  -198    53   108  -175  -225  -157  -186  -125   -54    82   151   -55  -125  -247  -213   -80    71   111
 I   -119    80  -379  -329   -39  -319  -221   250  -290   153    -2  -268  -307  -250  -282  -223  -123   121  -177  -151  -223  -330  -300  -138    71   111
 Y   -359  -274  -399  -424    60  -394   -88  -276  -369  -180  -227  -283  -396  -300  -343  -330  -343  -289    -1   473  -330  -349  -378  -278    71   111
 *   1364   426   528   746   801   595   353  1285   831  1837   573   686   535   608   750  1157  1110  1552   256   629     0     0     0     0 
*/
    
  
		
    List<String> lines  = FileUtils.readLines(new File(path));
    LinkedList<Double[]> input = new LinkedList();
	//G    -19  -134   -61    43  -171   223   -36  -144    -2  -140   -74     5  -166     6   -54     4   -26  -101   -42  -124     4   -32    29   -42   100   100

	
	for(int q= 11; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim();

		if(!one_line.startsWith("!")){

			String next_line = "";
			String[] split_str = space_pattern.split(one_line);

			if(split_str.length == 27){

				//
				Double[] arr = new Double[26];

				for(int w = 1; w < split_str.length; w++){

					arr[w-1] = new Double(split_str[w].trim());

				}

				
				input.add(arr);
			}//20
		}
		
			
 	}//q

	System.out.println("line size:"+ pfam_id + ":" + input.size() );

	return(input);
 

}//method


static LinkedHashMap<String,LinkedList<String>> pdbid_to_pfam;
static LinkedHashMap<LinkedList<String>,LinkedList<String>> pfam_to_pdbid;

void map_pdbid_to_pfam() throws IOException{

	System.out.println("map_pdbid_to_pfam()" );

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");

	pdbid_to_pfam = new LinkedHashMap();
	pfam_to_pdbid = new LinkedHashMap();

	//the number of the compound type was 1881

	String path = "data/sql.pdbid.pfam.final";

	List<String> lines  = FileUtils.readLines(new File(path));


	for(int q= 0; q < lines.size(); q++){

		/*
 101m Globin
102l Phage_lysozyme
102m Globin
6zym WD40_LIG_16-1
7dvq WD40_LIG_16-0
7dvq WD40_LIG_16-1
*/

  		String one_line = lines.get(q).trim();
		
		String[] tab_split = space_pattern.split(one_line);

		String pdb_name = tab_split[0].trim();
		String domain_name = tab_split[1].trim();
		String pfam = domain_to_pfam.get(domain_name);
	
		if(pdbid_to_pfam.containsKey(pdb_name)){
			
				LinkedList<String> exist = pdbid_to_pfam.get(pdb_name);
				exist.add(pfam);
				pdbid_to_pfam.put(pdb_name,exist);
		}else{
			
				LinkedList<String> exist = new LinkedList();
				exist.add(pfam);
				pdbid_to_pfam.put(pdb_name,exist);
		}
	}// q lines

	List<String> keys_list = new LinkedList<String>(pdbid_to_pfam.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		String pdb = keys_list.get(i);
		LinkedList<String> pfam = pdbid_to_pfam.get(pdb);
		//System.out.println("error:" + pdb + ":" );
		//System.out.println("error:" + Arrays.toString(pfam.toArray()) + ":" + pfam.size());
		
		boolean null_check=false;
		
		for(int y = 0; y < pfam.size(); y++){
		
			if(pfam.get(y) == null){
			
				null_check =true;
				break;
			}
		}
		
		if(pfam != null && !Arrays.toString(pfam.toArray()).contains("null") && !null_check){
		
			//System.out.println("error inside:" + Arrays.toString(pfam.toArray()) + ":" );
			Collections.sort((List)pfam);

	//if(comp.containsAll(conden_cat_comp) && conden_cat_comp.containsAll(comp) && comp.size() == conden_cat_comp.size()){

			if(pfam_to_pdbid.containsKey(pfam)){
				
					LinkedList<String> exist = pfam_to_pdbid.get(pfam);

					if(!exist.contains(pdb)){
						exist.add(pdb);
					}
					pfam_to_pdbid.put(pfam,exist);
			}else{
				
					LinkedList<String> exist = new LinkedList();
					exist.add(pdb);
					pfam_to_pdbid.put(pfam,exist);
			}
		}
	}

	
	
}//method

static LinkedHashMap<String,LinkedList<String>> pdb_to_uniprot;
static LinkedHashMap<String,LinkedList<String>> uniprot_to_pdb;

void  map_pdb_to_uniprot()  throws IOException{

	pdb_to_uniprot = new LinkedHashMap();
	uniprot_to_pdb = new LinkedHashMap();

	String path = "data/uniprot-pdb.tsv.only.uni.pdb";
 	
	pdb_to_uniprot = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern left_pattern = Pattern.compile("left");
	Pattern semi_pattern = Pattern.compile(";");

	List<String> lines  = FileUtils.readLines(new File(path));
	
	//A0A022MRT4 6SIW;6SIX;6SIY;6SIZ;6TM4;				
	for(int k = 0; k <lines.size(); k++){
				
		String one_line = lines.get(k).toUpperCase();
		String[] space_split = space_pattern.split(one_line);
		String uniprot = space_split[0].trim();
		String pdbs = space_split[1].trim();
		String[] semi_split = semi_pattern.split(pdbs);
		LinkedList<String> pdb_list = new LinkedList();

		for(int q = 0; q < semi_split.length; q++){

			String pdb = semi_split[q].trim();
			pdb_list.add(pdb);
			if(pdb_to_uniprot.containsKey(pdb)){

				LinkedList<String> exist = pdb_to_uniprot.get(pdb);
				if(!exist.contains(uniprot)){
					exist.add(uniprot);
				}
				pdb_to_uniprot.put(pdb,exist);

			}else{
				LinkedList<String> exist = new LinkedList();
				exist.add(uniprot);
				
				pdb_to_uniprot.put(pdb,exist);
			}
		}
		uniprot_to_pdb.put(uniprot,pdb_list);
	}//k

	System.out.println("uniprot_to_pdb:"+uniprot_to_pdb.size());
	System.out.println("pdb_to_uniprot:"+pdb_to_uniprot.size());
		
    }

static LinkedHashMap<String,LinkedList<String>> atxg_to_pfam;
static LinkedHashMap<LinkedList<String>,LinkedList<String>> nr_pfam_to_atxg;
static LinkedList<LinkedList<String>> pfam_tair_gene_list;
 
void map_pfam_to_tair_gene() throws IOException{

 	System.out.println("map_pfam_to_tair_gene()" );
 
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	
	atxg_to_pfam= new LinkedHashMap();
	nr_pfam_to_atxg= new LinkedHashMap();
	pfam_tair_gene_list= new LinkedList();
	
 /*
 	String path = "/run/media/ubuntu2/data052023/art_foundation_log/art_foundation/db/out_get_pfam_from_genbank_tair";
	List<String> lines  = FileUtils.readLines(new File(path));
					
	for(int k = 0; k < lines.size(); k++){

				
		String one_line = lines.get(k).toUpperCase();
		String[] colon_split = colon_pattern.split(one_line);
		//System.out.println("length check:" + colon_split.length);
		String atxg = colon_split[0].trim();
		String pfam_str = colon_split[colon_split.length-1].trim();
		LinkedList<String> pfam_list = new LinkedList();
		
		if(pfam_str.contains(",")){
		
			String[] comma_split = comma_pattern.split(pfam_str);
			for(int z = 0; z < comma_split.length; z++){
				pfam_list.add(comma_split[z].trim());
			}
				
		}else{
		
			pfam_list.add(pfam_str);
		}
		Collections.sort((List)pfam_list);
		
		if(atxg_to_pfam.containsKey(atxg)){
		
			LinkedList<String> exist = atxg_to_pfam.get(atxg);
			exist.addAll(pfam_list);
			
			Set<String> pfams_h_inv = new LinkedHashSet<>();
			pfams_h_inv.addAll((List)exist);
					
			// Clear the list
			((List)exist).clear();
											  
			// add the elements of set
			// with no duplicates to the list
			((List)exist).addAll(pfams_h_inv);
			Collections.sort((List)exist);
			
			
			atxg_to_pfam.put(atxg,exist);
			
		}else{
		
			atxg_to_pfam.put(atxg,pfam_list);
		}
				
		
	}//k
	*/

	
	
 
 	String path = "data/out_mix_interpro_genbank.final";
	//AT1G01030:[PF02362]
	//AT1G01040:[PF00035, PF00271, PF00636, PF02170, PF04851, PF14622]
	
	
	List<String> lines  = FileUtils.readLines(new File(path));
	
					
	for(int k = 0; k <lines.size(); k++){
				
		String one_line = lines.get(k);

		String[] colon_split = colon_pattern.split(one_line);
		String atxgid = colon_split[0].trim();
		String pfamid = colon_split[1].trim();
		pfamid=pfamid.replaceAll("\\[","");
		pfamid=pfamid.replaceAll("\\]","");
		//System.out.println("atxgid:"+atxgid);
		//System.out.println("pfamid:"+pfamid);

		LinkedList<String> nr_pfam = new LinkedList();
		String[] comma_split = comma_pattern.split(pfamid);
		for(int i = 0; i < comma_split.length; i++){

			nr_pfam.add(comma_split[i].trim());
		}

		atxg_to_pfam.put(atxgid,nr_pfam);

		//redundant
			
		Set<String> pfams_h_inv = new LinkedHashSet<>();
		pfams_h_inv.addAll((List)nr_pfam);
									
		// Clear the list
		((List)nr_pfam).clear();
															  
		// add the elements of set
		// with no duplicates to the list
		((List)nr_pfam).addAll(pfams_h_inv);
		Collections.sort((List)nr_pfam);
					

		if(nr_pfam_to_atxg.containsKey(nr_pfam)){
			
			LinkedList<String> exist = nr_pfam_to_atxg.get(nr_pfam);
				
			if(!exist.contains(atxgid)){
				
				exist.add(atxgid);
			}
				
			nr_pfam_to_atxg.put(nr_pfam, exist);
			
		}else{
			
			LinkedList<String> exist = new LinkedList();
				
			exist.add(atxgid);
			nr_pfam_to_atxg.put(nr_pfam, exist);
		}
				
	}//k
		
	
	
	
	List<String> keys_list = new LinkedList<String>(atxg_to_pfam.keySet());
	
	for(int i = 0; i < keys_list.size(); i++){
	
		pfam_tair_gene_list.add(atxg_to_pfam.get(keys_list.get(i)));
	}
	
	System.out.println("atxg_to_pfam.size():" + atxg_to_pfam.size());

	

}//methodmap_pfam_human_gene

static LinkedHashMap<String,LinkedList<String>> tair_gene_uniprot;
static LinkedHashMap<String,LinkedList<String>> uniprot_tair_gene;
	

void map_atxg_to_uniprot() throws IOException{

	System.out.println("map_atxg_to_uniprot():" );

	tair_gene_uniprot = new LinkedHashMap();
	uniprot_tair_gene = new LinkedHashMap();

  	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	

	String input="data/TAIR2UniprotMapping-JAN2023.txt.col1.2";
	
	List<String> lines  = FileUtils.readLines(new File(input));

	
	for(int q= 0; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim().toUpperCase();
		
		if(!one_line.equals("")){

			

			String[] tab_split = space_pattern.split(one_line);
			String uniprot = tab_split[0].trim();
	
			String tairid = tab_split[1].trim();

			if(tair_gene_uniprot.containsKey(tairid)){
				LinkedList<String> exist = tair_gene_uniprot.get(tairid);
				exist.add(uniprot);
				tair_gene_uniprot.put(tairid,exist);

			}else{
				LinkedList<String> exist = new LinkedList();
				exist.add(uniprot);
				tair_gene_uniprot.put(tairid,exist);

			}

			if(uniprot_tair_gene.containsKey(uniprot)){
				LinkedList<String> exist = uniprot_tair_gene.get(uniprot);
				exist.add(tairid);
				uniprot_tair_gene.put(uniprot,exist);

			}else{
				LinkedList<String> exist = new LinkedList();
				exist.add(tairid);
				uniprot_tair_gene.put(uniprot,exist);

			}
			
	
			LinkedList<String> array = new LinkedList();
			array.add(uniprot);
			tair_gene_uniprot.put(tairid,array);


		}
	
	}//q
	System.out.println("tair_gene_uniprot:" + tair_gene_uniprot.size());
	
}

static LinkedHashMap<String,LinkedList<String>> string_db_tair;

void make_string_db_tair() throws IOException{
	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	string_db_tair = new LinkedHashMap();

	
	String path = "big/3702.protein.physical.links.detailed.v11.5.txt.no3702.cutoff.500.atxg.only";

    	
		
    List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){
	
	
		String one_line = lines.get(q).trim();
		String[] split_str = colon_pattern.split(one_line);

		String tf_gene = split_str[0];
		String cofactor_gene = split_str[1];

		if(string_db_tair.containsKey(tf_gene)){
					
			LinkedList<String> exist = string_db_tair.get(tf_gene);
			if(!exist.contains(cofactor_gene)){
							
				exist.add(cofactor_gene);
			}
			string_db_tair.put(tf_gene,exist);
												
		}else{
					
			LinkedList<String> exist = new LinkedList();
			exist.add(cofactor_gene);
			string_db_tair.put(tf_gene,exist);
						
		}
					
			
		
			
 	}//q
 	
 }	

static LinkedHashMap<String,LinkedList<String>> string_db_tf_co_tair;

void make_tf_cofactor_list_tair() throws IOException{
	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	string_db_tf_co_tair = new LinkedHashMap();

	
	String path = "big/3702.protein.physical.links.detailed.v11.5.txt.no3702.cutoff.500.atxg.only";

    	List<String> tf_lines  = FileUtils.readLines(new File("data/TF.atgx.sort.uniq"));
    
    
  
		
    List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){
	
	
		String one_line = lines.get(q).trim();
		String[] split_str = colon_pattern.split(one_line);

		String tf_gene = split_str[0];
		String cofactor_gene = split_str[1];
		
		
		if(tf_lines.contains(tf_gene)){
				
					if(string_db_tf_co_tair.containsKey(tf_gene)){
					
						LinkedList<String> exist = string_db_tf_co_tair.get(tf_gene);
						if(!exist.contains(cofactor_gene)){
						
							exist.add(cofactor_gene);
						}
						string_db_tf_co_tair.put(tf_gene,exist);
												
					}else{
					
						LinkedList<String> exist = new LinkedList();
						exist.add(cofactor_gene);
						string_db_tf_co_tair.put(tf_gene,exist);
						
					}
					
			
		}else if(tf_lines.contains(cofactor_gene)){
		
				
						
					if(string_db_tf_co_tair.containsKey(cofactor_gene)){
						
							LinkedList<String> exist = string_db_tf_co_tair.get(cofactor_gene);
							if(!exist.contains(tf_gene)){
							
								exist.add(tf_gene);
							}
							string_db_tf_co_tair.put(cofactor_gene,exist);
							
							
					}else{
						
							LinkedList<String> exist = new LinkedList();
							exist.add(tf_gene);
							string_db_tf_co_tair.put(cofactor_gene,exist);
							
							
					}
						
		}
			
 	}//q
 	
 }	


}//inner class ddi



static class pathway_anal_module{



static String[] species = {"arabidopsis_thaliana"};
//,"brachypodium_distachyon","citrullus_lanatus","cucumis_sativus","gossypium_raimondii","malus_domestica","oryza_sativa","phaseolus_vulgaris","populus_trichocarpa","prunus_persica","setaria_italica","solanum_lycopersicum","solanum_tuberosum","sorghum_bicolor","vitis_vinifera","zea_mays"};

static String[] spec_genome = {"Arabidopsis@TAIR10"};
//,"Brachypodium@v3.0","Citrullus@v1.0","Cucumis@ASM407-v2","Gossypium@v2.1","Malus@GDDH13 v1.1","Oryza@IRGSP-1.0","Phaseolus@v2.1","Populus@v3.0","Prunus@v2.0","Setaria@v2.0","Solanum@ITAG2.4","Solanum@v4.03","Sorghum@v3.0.1","Vitis@IGGP_12X","Zea@AGPv3"};

static String[] gene_name_format = {"AT1G01030"};
//,"Bradi1g00200","Cla000003.g","Csa_1G000010","Gorai.001G000100","gene:MD00G1000600","LOC_Os01g01010","Phvul.001G000700","Potri.001G000800","Prupe.1G000400","SETIT_018185mg","Solyc00g005040.2","PGSC0003DMG400039401","Sobic.001G001800","GSVIVG01012259001","AC177838.2_FG015"};
    
static LinkedHashMap<String,String> map_compound_type_to_type_group;
    
//retrieve compound type group
void retrieve_plantcyc_comp_type() throws IOException{

	//System.out.println("retrieve_plantcyc_comp_type:" );
	map_compound_type_to_type_group = new LinkedHashMap();

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");

	String comp_all_path = "data/plantcyc/plantcyc/comp.types.sort.uniq";

	String first_path = "data/plantcyc/plantcyc/compound.type.check.first";

	String comp_type_path = "data/plantcyc/plantcyc/compound_type";

	String last_check_path = "data/plantcyc/plantcyc/last_check";

	String last_first_path = "data/plantcyc/plantcyc/last_first";

	List<String> all_lines  = FileUtils.readLines(new File(comp_all_path));
	List<String> first_lines  = FileUtils.readLines(new File(first_path));
	List<String> comp_type_lines  = FileUtils.readLines(new File(comp_type_path));
	List<String> last_check_lines  = FileUtils.readLines(new File(last_check_path));
	List<String> last_first_lines  = FileUtils.readLines(new File(last_first_path));

	LinkedList<String> remain = new LinkedList();
	LinkedList<String> first = new LinkedList();
	LinkedList<String> type = new LinkedList();
	LinkedList<String> last_check = new LinkedList();
	LinkedList<String> last_first = new LinkedList();

	for(int q= 0; q < first_lines.size(); q++){

  		String one_line = first_lines.get(q).trim().toUpperCase();
		first.add(one_line);
	}

	for(int q= 0; q < comp_type_lines.size(); q++){

  		String one_line = comp_type_lines.get(q).trim().toUpperCase();
		type.add(one_line);
	}
	for(int q= 0; q < last_check_lines.size(); q++){

  		String one_line = last_check_lines.get(q).trim().toUpperCase();
		last_check.add(one_line);
	}
	for(int q= 0; q < last_first_lines.size(); q++){

  		String one_line = last_first_lines.get(q).trim().toUpperCase();
		last_first.add(one_line);
	}

	for(int q= 0; q < all_lines.size(); q++){

		String category = "";

  		String one_line = all_lines.get(q).trim().toUpperCase();
  		one_line = one_line.replaceFirst("TYPES - ","").trim();
		boolean first_flag = false;

		//System.out.println("498 one_line:" + one_line);

		for(int r = 0; r < first.size(); r++){

			if(one_line.contains(first.get(r))){

				first_flag = true;
				category = "F_" + first.get(r);
				//System.out.println(" :" + one_line + ":category f:" + category);
				break;
			}
		}

		//System.out.println("first_flag:" + first_flag);

		if(!first_flag){

			String[] dash_split = dash_pattern.split(one_line);

			//System.out.println("dash_split[dash_split.length-1]:" + dash_split[dash_split.length-1]);

			String last_dash = dash_split[dash_split.length-1];
			if(last_dash.substring(last_dash.length()-1, last_dash.length()).equals("S")){

				last_dash =last_dash.substring(0,last_dash.length()-1);
			}

			if(type.contains(last_dash)){

	
				for(int a = 0; a < type.size(); a++){

					if(last_dash.equals(type.get(a))){
						category = "D_" + type.get(a);
						//System.out.println("category d:" + category);
						break;
						
					}
				}//a

			}else{
				//System.out.println("!type.contains:" + dash_split[dash_split.length-1]);

				boolean check = false;

				for(int e = 0; e < type.size(); e++){

					if(one_line.contains(type.get(e))){

						//System.out.println("type.get(e):" + type.get(e));

						category = "A_" + type.get(e);
						
						check = true;
						break;

					}
				}//e

				//System.out.println("check:" + check);


				if(!check){

					boolean last_flag = false;
					for(int w = 0; w < last_first.size(); w++){

						if(one_line.contains(last_first.get(w))){
							category = "T_" + last_first.get(w);
							//System.out.println("category t:" + category);
							last_flag = true;
							break;
						}
					}

					if(!last_flag){

						for(int w = 0; w < last_check.size(); w++){

							String temp = last_check.get(w).trim();
							int temp_length = temp.length();
							if(one_line.substring(one_line.length()-1,one_line.length()).equals("S")){

								one_line = one_line.substring(0, one_line.length()-1);
								//System.out.println("category s:" + one_line);
							}
							
							if(one_line.length()-temp_length >0 && one_line.length()-temp_length <one_line.length()){
								String last_word =  one_line.substring(one_line.length()-temp_length,one_line.length() );

								if(temp.equals(last_word)){

									category = "L_" + last_check.get(w);
									//System.out.println("category l:" + category);
									break;
								}
							}
						}//w
					}//flag
				}//!check
			}//else
		}//if not first
		
		if(category.equals("")){
			category = "NA";
		}
		
		map_compound_type_to_type_group.put(one_line,category);
		//System.out.println("type group final:" + one_line + ":" + category);
	}// q lines
	
}//method


//pmn compounds.dat file

static String[] struc_group =  {"A-DARA-HEX-2:5_2:KETO","A-DARA-PEN-1:4","A-DGAL-HEX-1:4","A-DGAL-HEX-1:5","A-DGAL-HEX-1:5_6:A","A-DGLC-HEX-1:4_6:A","A-DGLC-HEX-1:5","A-DGLC-HEX-1:5_6:A","A-DGLC-HEX-1:5_6:D","A-DGRO-DGAL-NON-2:6_1:A_2:KETO_3:D","A-DIDO-HEX-1:5_6:A","A-DMAN-HEX-1:5","A-DMAN-HEX-1:5_6:D","A-DMAN-OCT-2:6_1:A_2:KETO_3:D","ADP","A-DXYL-PEN-1:5","A-HEP-1:5","A-HEX-1:5","A-HEX-1:5_6:A","A-LARA-PEN-1:4","A-LGAL-HEX-1:5","A-LGAL-HEX-1:5_6:D","A-LIDO-HEX-1:5_6:A","A-LMAN-HEX-1:5_6:D","A-LXYL-PEN-1:5","A-XGAL-HEX-1:5_6:A","B-DARA-HEX-2:5_2:KETO","B-DGAL-HEX-1:4","B-DGAL-HEX-1:5","B-DGAL-HEX-1:5_6:A","B-DGLC-HEX-1:5","B-DGLC-HEX-1:5_6:A","B-DGLC-HEX-1:5_6:D","B-DGRO-DGAL-NON-2:6_1:A_2:KETO_3:D","B-DIDO-HEX-1:5_6:A","B-DMAN-HEX-1:5","B-DMAN-HEX-1:5_6:A","B-DMAN-HEX-1:5_6:D","B-DMAN-OCT-2:6_1:A_2:KETO_3:D","B-DXYL-PEN-1:5","B-HEX-1:5","B-HEX-1:5_6:A","B-LARA-PEN-1:4","B-LGAL-HEX-1:5_6:D","B-LGLC-HEX-1:5_6:A","B-LIDO-HEX-1:5_6:A","B-LMAN-HEX-1:5_6:D","B-PEN-1:5","B-XXYL-PEN-1:5","Ceramides","CMP","CPD0-1283","CPD-23878","DOLICHOLP","EGF-Domain-Protein-L-Serine","FERULIC-ACID","GDP","GLC-D-LACTONE","GLUCONATE","GLYCO-CT-ACETYL","GLYCO-CT-AMINO","GLYCO-CT-METHYL","GLYCO-CT-N-ACETYL","GLYCO-CT-N-SULFATE","GLYCO-CT-PHOSPHATE","GLYCO-CT-PYR","GLYCO-CT-SULFATE","Lipid-A","LIPID-IV-A","Protein-hydroxyprolines","Protein-L-Asparagine","Protein-L-serine-or-L-threonine","Protein-L-serines","R","R1","R2","R3","TDP","THR","UDP","UNDECAPRENYL-DIPHOSPHATE","X-DARA-HEX-2:5_2:KETO","X-DARA-PEN-1:4","X-DGAL-HEX-1:5","X-DGLC-HEX-1:5","X-DGLC-HEX-1:5_6:A","X-DGLC-HEX-1:5_6:D","X-DMAN-HEX-1:5","X-DMAN-HEX-X:X","X-DXYL-PEN-1:5","X-HEP-1:5","X-HEP-X:X","X-HEX-1:5","X-LARA-PEN-1:4","X-LGAL-HEX-1:5_6:D","X-LIDO-HEX-1:5_6:A","X-LMAN-HEX-1:5_6:D","X-PEN-1:4"};

static String[] formula_ele = {"AG","AL","AS","B","BA","BE","BR","C","CA","CD","CL","CO","CS","CU","F","FE","H","HG","I","K","LI","MG","MN","MO","N","NA","NI","O","P","PB","PT","RB","S","SE","SR","TE","V","W","ZN"};

static String[] comp_type_group = {"A_ADENINE","A_ALCOHOL","A_ALDEHYDE","A_ALDODIOSE","A_ALDOSE","A_ALKALOID","A_ALKANE","A_AMINE","A_AMINO-ACID","A_ANTIGEN","A_BENZOATE","A_BROMIDE","A_CAMPHOR","A_CARBOXYLATE","A_CHLORIDE","A_CHOLINE","A_CYANIDIN","A_CYTIDINE","A_CYTOSINE","A_DIOL","A_EPOME","A_EPOXIDE","A_ESTER","A_ETHER","A_FATTY-ACID","A_FLAVANONE","A_FLAVIN","A_FLAVONE","A_FLAVONOID","A_FOLATE","A_FRUCTOSE","A_FURANOSE","A_GANGLIOSIDE","A_GLUCAN","A_GLUCOSE","A_GLUTATHIONE","A_GLYCAN","A_GLYCERATE","A_GLYCERIDE","A_GLYCEROL","A_GLYCO-","A_HEME","A_HEXOSE","A_HORMONE","A_IMINE","A_INDOLE","A_KETONE","A_LAMINE","A_LIPID","A_NITRILE","A_N-TERMINAL","A_NUCLEOTIDES","A_OXIME","A_PEPTIDE","A_PHENOL","A_PHOSPHATE","A_POLYMER","A_PURINE-","A_PYRIMIDINE-","A_QUINOL","A_QUINONE","A_RNA-","A_SACCHARIDE","A_STEROID","A_STEROL","A_SUGAR","A_SULFONATE","A_SULFUR","A_UDP-","A_URACIL","A_UREA","A_URIDINE","A_VITAMIN","A_XYLAN","D_ADENINE","D_ALCOHOL","D_ALDEHYDE","D_ALDOSE","D_ALKALOID","D_ALKANE","D_AMINE","D_ANTIGEN","D_CARBOXYLATE","D_COA","D_CYTIDINE","D_CYTOSINE","D_DIOL","D_EPOXIDE","D_ESTER","D_ETHER","D_FLAVANONE","D_FLAVIN","D_FLAVONE","D_FLAVONOID","D_FOLATE","D_FRUCTOSE","D_GANGLIOSIDE","D_GLUCAN","D_GLUCOSE","D_GLUTATHIONE","D_GLYCAN","D_GLYCEROL","D_HEXOSE","D_HORMONE","D_KETONE","D_LIPID","D_NITRILE","D_PEPTIDE","D_PHENOL","D_PHOSPHATE","D_POLYMER","D_QUINOL","D_QUINONE","D_STEROID","D_STEROL","D_SUGAR","D_SULFATE","D_SULFONATE","D_URACIL","D_URIDINE","D_VITAMIN","D_XYLAN","F_ABSCISIC-ACID","F_AUXIN","F_COFACTOR","F_CYTOKININ","F_DNA","F_GIBBERELLINS","F_JASMONATE","F_LIGAND","F_NAD","F_NDP","F_PROTEIN","F_RRNA","F_TRNA","L_ACID","L_AL","L_ALKENE","L_AM","L_AN","L_ATE","L_EM","L_EN","L_ENAL","L_ENE","L_ETE","L_IN","L_INE","L_ION","L_MIA","L_NAL","L_NONE","L_NOSE","L_ODE","L_OID","L_OL","L_OLE","L_OME","L_ONE","L_OR","L_ORE","L_OSE","L_OXIDE","L_PYRAN","L_RIN","L_RING","L_SAN","L_THIOL","L_TONE","L_TOSE","L_TRIOSE","L_YLL","T_ACCEPTOR","T_ANION","T_ANTIBIOTIC","T_CARRIER","T_CATION","T_COMPLEXE","T_COMPOUND","T_CPD","T_DERIVATIVE","T_ENZYME","T_FACTOR"};

/*
UNIQUE-ID - 1-PHOSPHATIDYL-1D-MYO-INOSITOL-35-BISPH
TYPES - Phosphoinositides
COMMON-NAME - 1-phosphatidyl-1D-myo-inositol 3,5-bisphosphate
ATOM-CHARGES - (24 -1)
ATOM-CHARGES - (22 -1)
ATOM-CHARGES - (21 -1)
ATOM-CHARGES - (19 -1)
ATOM-CHARGES - (18 -1)
COMMENT - The phosphorylated derivatives of |FRAME: L-1-phosphatidyl-inositols "phosphatidylinositol"| (PtdIns), collectively called phosphoinositides (PIs), are eukaryotic cell membrane phospholipids whose inositol head group is phosphorylated at position D3, D4, or D5, singly or in all combinations, to yield seven PI metabolites.

*/

//UNIQUE-ID
//TYPES
//CHEMICAL-FORMULA
//MOLECULAR-WEIGHT
//STRUCTURE-GROUPS

static LinkedHashMap<String,String>[] compound_id_to_struc_group;
static LinkedHashMap<String,Boolean[]>[] compound_id_to_formula_ele;
static LinkedHashMap<String,Double>[] compound_id_to_mw;
static LinkedHashMap<String,String>[] compound_id_to_type;
static LinkedHashMap<String,String>[] compound_id_to_type_group;

void map_compound_id_to_features() throws IOException{

	//System.out.println("map_compound_id_to_features():" );

	compound_id_to_struc_group = new LinkedHashMap[species.length];
	compound_id_to_formula_ele = new LinkedHashMap[species.length];
	compound_id_to_mw = new LinkedHashMap[species.length];
	compound_id_to_type = new LinkedHashMap[species.length];
	compound_id_to_type_group = new LinkedHashMap[species.length];
	
	String[] species_file_path = {"big/pmn/aracyc/17.1.2/data/"};
	//,"big/pmn/brachypodiumcyc/8.0.0/data/","big/pmn/watermeloncyc/1.0.2/data/","big/pmn/csativa_pkcyc/3.0.2/data/","big/pmn/graimondiicyc/3.0.2/data/","big/pmn/mdomesticacyc/3.0.2/data/","big/pmn/oryzacyc/7.2.1/data/","big/pmn/commonbeancyc/3.0.2/data/","big/pmn/poplarcyc/13.0.0/data/","big/pmn/ppersicacyc/3.0.2/data/","big/pmn/setariacyc/7.1.1/data/","big/pmn/tomatocyc/6.0.0/data/","big/pmn/potatocyc/7.0.0/data/","big/pmn/sorghumbicolorcyc/7.1.1/data/","big/pmn/grapecyc/10.0.0/data/","big/pmn/corncyc/11.0.0/data/"};

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");

	for(int i = 0; i < species_file_path.length; i++){

		String path = species_file_path[i] + "compounds.dat";
		compound_id_to_struc_group[i] = new LinkedHashMap();
		compound_id_to_formula_ele[i] = new LinkedHashMap();
		compound_id_to_mw[i] = new LinkedHashMap();
		compound_id_to_type[i] = new LinkedHashMap();
		compound_id_to_type_group[i] = new LinkedHashMap();
	
	
		List<String> lines  = FileUtils.readLines(new File(path));
		
		/*
		
		UNIQUE-ID - Phosphatidylserines-34-3
TYPES - L-1-PHOSPHATIDYL-SERINE
COMMON-NAME - phosphatidylserine-34:3
COMMENT - In this group of Phosphatidylserines, the length of the two acyl chains, and the location and number of the double bonds varies. For each of the group members, the total number of carbon atoms in the two acyl chains adds up to 34, and the total number of double bonds in the two acyl groups adds up to 3.
HAS-NO-STRUCTURE? - T
INSTANCE-NAME-TEMPLATE - CPD-*
SYNONYMS - PS-34:3
SYNONYMS - PS34:3


		*/
		
		Boolean[] form_temp = new Boolean[formula_ele.length];
		for(int q= 0; q<formula_ele.length; q++){
		
			form_temp[q] = new Boolean(false);
		
		}

		int count = 0;
		String id = "";
		String type = "";
		String type_group = "";
		String chem_form = "";
		String mw = "";
		String struc_group = "";

		for(int q= 0; q < lines.size(); q++){
		
			String one_line = lines.get(q).trim().toUpperCase().trim();
			
			
			
			
			if(!one_line.equals("") && !one_line.startsWith("#")){
			
			
				if(one_line.contains("UNIQUE-ID ")){
				
					if(count != 0){
					
						compound_id_to_formula_ele[i].put(id,form_temp);
						//System.out.println("error 706:" + id + ":" + Arrays.toString(form_temp));
						form_temp = new Boolean[formula_ele.length];
						for(int e= 0;e< form_temp.length; e++){
						
							form_temp[e] = new Boolean(false);
						
						}
						
						id = "";
						type = "";
						type_group = "";
						chem_form = "";
						mw = "";
						struc_group = "";
					}
				
					String[] space_split = space_pattern.split(one_line);
					id = space_split[space_split.length-1].trim();
					count = count +1;
				
				}else if(one_line.contains("TYPES ")){
				
					String[] space_split = space_pattern.split(one_line);
					type = space_split[space_split.length-1].trim();
					type_group = map_compound_type_to_type_group.get(type);
					//System.out.println("724:" + id + ":" + type + ":" + type_group);
					compound_id_to_type[i].put(id,type);
					compound_id_to_type_group[i].put(id,type_group);
				
				
				}else if(one_line.contains("CHEMICAL-FORMULA ")){
				
					String comp_temp = one_line.replaceFirst("CHEMICAL-FORMULA - ","");
				
					for(int z = 0; z < form_temp.length; z++){
					
						if(comp_temp.contains(formula_ele[z])){
						
							form_temp[z] = new Boolean(true);
							//System.out.println("737:" + id + ":true" );
							break;
						
						}
					}
		
				
				}else if(one_line.contains("MOLECULAR-WEIGHT ")){

					//System.out.println("mw:"+ one_line);
				
					String[] space_split = space_pattern.split(one_line);

					if(space_split.length == 3){
						mw = space_split[space_split.length-1].trim();
						compound_id_to_mw[i].put(id,new Double(mw));
						//System.out.println("753:"+ id + ":" + mw);
					}else{
						compound_id_to_mw[i].put(id,Double.NaN);
						//System.out.println("756:"+ id + ":Double.NaN" );

					}
					
					
				
				
				}else if(one_line.contains("STRUCTURE-GROUPS ")){
				
					String[] space_split = space_pattern.split(one_line);
					struc_group = space_split[space_split.length-1].trim();
					//System.out.println("struc_group:"+ struc_group);
					String[] dash_split = dash_pattern.split(struc_group);
					String group_temp = "";
					for(int z = 0; z < dash_split.length; z++){
					
						if(dash_split[z].trim().matches(".*\\d.*")){
						
							break;
						
						}
						
						group_temp = group_temp + dash_split[z].trim() + "-";
						
						//System.out.println("779:"+ group_temp);
					}

					if(group_temp.equals("")){
						group_temp = "NA";
						//System.out.println("784:"+ group_temp);
					}else{
						
						group_temp = group_temp.substring(0, group_temp.length()-1);
						//System.out.println("789:"+ group_temp);
					}
					//System.out.println("791:"+ id + ":" + group_temp);
					compound_id_to_struc_group[i].put(id,group_temp);
				
				}

				
				
			}
		
		}//q

	}

}



//reactions
static String[] rxn_type = {"Chemical-Reactions","Complex-Processes","Composite-Reactions","Electron-Transfer-Reactions","Light-Interactions","Polysaccharide-Reactions","Protein-Modification-Reactions","Protein-Reactions","Redox-Half-Reactions","RNA-Reactions","Small-Molecule-Reactions","TR-11","TR-12","TR-13","Transport-Reactions","tRNA-Charging-Reactions","tRNA-Reactions","Unknown-Conversions"};

static String[] pwy_list = {"DETOX1-PWY-1","PWY-1001","PWY-101","PWY-102","PWY-1042","PWY-1061","PWY-1081","PWY-1121","PWY-116","PWY-1269","PWY-1422","PWY-1581","PWY-1801","PWY-181","PWY-1822","PWY-1881","PWY-2","PWY-2161","PWY-2161B","PWY-2161B-PMN","PWY-2261","PWY-2301","PWY-2501","PWY-2541","PWY-2582","PWY-2681","PWY-2781","PWY-282","PWY-2841","PWY-2881","PWY-2901","PWY-2902","PWY-2981","PWY-3101","PWY-3181","PWY-321","PWY-3221","PWY-3261","PWY-3341","PWY-3385","PWY-3461","PWY-3462","PWY-3542","PWY-3561","PWY-361","PWY-3742","PWY-3781","PWY-3801","PWY-381","PWY-3821","PWY-3841","PWY-3881","PWY-3941","PWY-3981","PWY-3982","PWY-4","PWY-401","PWY-4021","PWY-4041","PWY-4081","PWY-4101","PWY-4203","PWY-4261","PWY-43","PWY-4302","PWY-4321","PWY-4341","PWY-4361","PWY-4381","PWY-4421","PWY-4541","PWY-4562","PWY-46","PWY-4601","PWY-4621","PWY-4661","PWY-4702","PWY-4821","PWY-4841","PWY-4861","PWY-4983","PWY-4984","PWY-5027","PWY-5032","PWY-5034","PWY-5035","PWY-5036","PWY-5041","PWY-5046","PWY-5048","PWY-5049","PWY-5060","PWY-5063","PWY-5064","PWY-5068","PWY-5070","PWY-5080","PWY-5083","PWY-5084","PWY-5086","PWY-5097","PWY-5098","PWY-5107","PWY-5111","PWY-5113","PWY-5115","PWY-5116","PWY-5120","PWY-5122","PWY-5123","PWY-5125","PWY-5129","PWY-5136","PWY-5137","PWY-5138","PWY-5139","PWY-5142","PWY-5143","PWY-5147","PWY-5148","PWY-5160","PWY-5172","PWY-5175","PWY-5176","PWY-5188","PWY-5194","PWY-5203","PWY-5268","PWY-5269","PWY-5271","PWY-5272","PWY-5285","PWY-5305","PWY-5307","PWY-5312","PWY-5313","PWY-5326","PWY-5337","PWY-5340","PWY-5343","PWY-5350","PWY-5365","PWY-5381","PWY-5386","PWY-5389","PWY-5406","PWY-5407","PWY-5408","PWY-5409","PWY-5410","PWY-5434","PWY-5441","PWY-5453","PWY-5466","PWY-5467","PWY-5469","PWY-5473","PWY-5481","PWY-5484","PWY-5486","PWY-561","PWY-5651","PWY-5653","PWY-5659","PWY-5667","PWY-5669","PWY-5670","PWY-5686","PWY-5687","PWY-5690","PWY-5691","PWY-5692","PWY-5697","PWY-5698","PWY-5704","PWY-5723","PWY-5733","PWY-5754","PWY-5774","PWY-5783","PWY-5784","PWY-5791","PWY-5800","PWY-5805","PWY-5807","PWY-5837","PWY-5859","PWY-5863","PWY-5871","PWY-5876","PWY-5884","PWY-5886","PWY-5905","PWY-5912","PWY-5921","PWY-5925","PWY-5934","PWY-5936","PWY-5943","PWY-5944","PWY-5945","PWY-5946","PWY-5947","PWY-5957","PWY-5968","PWY-5971","PWY-5973","PWY-5975","PWY-5980","PWY-5983","PWY-5989","PWY-5995","PWY-5997","PWY-6","PWY-6010","PWY-6012","PWY-6013","PWY-6014","PWY-6019","PWY-6029","PWY-6030","PWY-6039","PWY-6040","PWY-6051","PWY-6064","PWY-6066","PWY-6068","PWY-6069","PWY-6115","PWY-6118","PWY-6120","PWY-6121","PWY-6124","PWY-6147","PWY-6163","PWY-6164","PWY-6196","PWY-6199","PWY-621","PWY-6219","PWY-622","PWY-6220","PWY-6232","PWY-6233","PWY-6235","PWY-6239","PWY-6257","PWY-6275","PWY-6287","PWY-6291","PWY-6295","PWY-6297","PWY-63","PWY-6303","PWY-6305","PWY-6333","PWY-6348","PWY-6351","PWY-6352","PWY-6363","PWY-6367","PWY-6369","PWY-641","PWY-6424","PWY-6440","PWY-6441","PWY-6442","PWY-6443","PWY-6457","PWY-6465","PWY-6466","PWY-6468","PWY-6473","PWY-6475","PWY-6482","PWY-6483","PWY-6494","PWY-6498","PWY-6502","PWY-6515","PWY-6520","PWY-6527","PWY-6535","PWY-6540","PWY-6543","PWY-6544","PWY-6546","PWY-6549","PWY-6554","PWY-6556","PWY-6577","PWY-6596","PWY-6599","PWY-66","PWY-6602","PWY-6605","PWY-6606","PWY-6607","PWY-6613","PWY-6618","PWY-6619","PWY-6624","PWY-6631","PWY-6663","PWY-6665","PWY-6673","PWY-6707","PWY-6710","PWY-6717","PWY-6724","PWY-6730","PWY-6733","PWY-6745","PWY-6754","PWY-6762","PWY-6773","PWY-6786","PWY-6787","PWY-6799","PWY-6801","PWY-6802","PWY-6803","PWY-6804","PWY-6806","PWY-6809","PWY-6823","PWY-6825","PWY-6837","PWY-6890","PWY-6898","PWY-6902","PWY-6908","PWY-6909","PWY-6910","PWY-6920","PWY-6922","PWY-6927","PWY-6930","PWY-6933","PWY-6936","PWY-6938","PWY-6949","PWY-695","PWY-6950","PWY-6952","PWY-6959","PWY-6963","PWY-6964","PWY-699","PWY-6990","PWY-701","PWY-702","PWY-7029","PWY-7036","PWY-7039","PWY-7060","PWY-7067","PWY-7068","PWY-7071","PWY-7101","PWY-7120","PWY-7137","PWY-7140","PWY-7141","PWY-7145","PWY-7158","PWY-7170","PWY-7176","PWY-7177","PWY-7183","PWY-7184","PWY-7185","PWY-7186","PWY-7187","PWY-7188","PWY-7189","PWY-7193","PWY-7197","PWY-7199","PWY-7204","PWY-7205","PWY-7206","PWY-7212","PWY-7213","PWY-7214","PWY-7219","PWY-7221","PWY-7224","PWY-7226","PWY-7227","PWY-7232","PWY-7250","PWY-7252","PWY-7253","PWY-7270","PWY-7280","PWY-7325","PWY-7343","PWY-7344","PWY-7346","PWY-735","PWY-7356","PWY-7363","PWY-7388","PWY-7394","PWY-7416","PWY-7417","PWY-7436","PWY-7450","PWY-7452","PWY-7458","PWY-7463","PWY-7464","PWY-7468","PWY-7477","PWY-7478","PWY-7481","PWY-7484","PWY-7487","PWY-7489","PWY-7497","PWY-7546","PWY-7560","PWY-7590","PWY-7601","PWY-762","PWY-7634","PWY-7683","PWY-782","PWY-7856","PWY-7859","PWY-7861","PWY-7897","PWY-7909","PWY-7985","PWY-7995","PWY-8001","PWY-801","PWY-8011","PWY-8027","PWY-8028","PWY-82","PWY-84","PWY-842","PWY-861","PWY-862","PWY-882","PWY-922"};

static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_rxn_type;
static LinkedHashMap<String,String>[] enzrxn_id_to_pwy;
static LinkedHashMap<String,LinkedList<String>>[] enzrxn_id_to_substrate_compound;
static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_substrate_form_ele;
static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_substrate_struc_group;
static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_substrate_type_group;
static LinkedHashMap<String,LinkedList<String>>[] enzrxn_id_to_product_compound;
static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_product_form_ele;
static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_product_struc_group;
static LinkedHashMap<String,Boolean[]>[] enzrxn_id_to_product_type_group;

void map_enzrxn_id_to_features() throws IOException{

	//System.out.println("map_enzrxn_id_to_features():" );

	enzrxn_id_to_rxn_type = new LinkedHashMap[species.length];
	enzrxn_id_to_pwy = new LinkedHashMap[species.length];
	enzrxn_id_to_substrate_compound = new LinkedHashMap[species.length];
	enzrxn_id_to_substrate_form_ele = new LinkedHashMap[species.length];
	enzrxn_id_to_substrate_struc_group = new LinkedHashMap[species.length];
	enzrxn_id_to_substrate_type_group = new LinkedHashMap[species.length];
	enzrxn_id_to_product_compound = new LinkedHashMap[species.length];
	enzrxn_id_to_product_form_ele = new LinkedHashMap[species.length];
	enzrxn_id_to_product_struc_group = new LinkedHashMap[species.length];
	enzrxn_id_to_product_type_group = new LinkedHashMap[species.length];
	
	
	String[] species_file_path = {"big/pmn/aracyc/17.1.2/data/"};
	//,"big/pmn/brachypodiumcyc/8.0.0/data/","big/pmn/watermeloncyc/1.0.2/data/","big/pmn/csativa_pkcyc/3.0.2/data/","big/pmn/graimondiicyc/3.0.2/data/","big/pmn/mdomesticacyc/3.0.2/data/","big/pmn/oryzacyc/7.2.1/data/","big/pmn/commonbeancyc/3.0.2/data/","big/pmn/poplarcyc/13.0.0/data/","big/pmn/ppersicacyc/3.0.2/data/","big/pmn/setariacyc/7.1.1/data/","big/pmn/tomatocyc/6.0.0/data/","big/pmn/potatocyc/7.0.0/data/","big/pmn/sorghumbicolorcyc/7.1.1/data/","big/pmn/grapecyc/10.0.0/data/","big/pmn/corncyc/11.0.0/data/"};

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	for(int i = 0; i < species_file_path.length; i++){

		String path = species_file_path[i] + "reactions.dat";
		enzrxn_id_to_rxn_type[i] = new LinkedHashMap();
		enzrxn_id_to_pwy[i] = new LinkedHashMap();
		enzrxn_id_to_substrate_compound[i] = new LinkedHashMap();
		enzrxn_id_to_substrate_form_ele[i] = new LinkedHashMap();
		enzrxn_id_to_substrate_struc_group[i] = new LinkedHashMap();
		enzrxn_id_to_substrate_type_group[i] = new LinkedHashMap();
		enzrxn_id_to_product_compound[i] = new LinkedHashMap();
		enzrxn_id_to_product_form_ele[i] = new LinkedHashMap();
		enzrxn_id_to_product_struc_group[i] = new LinkedHashMap();
		enzrxn_id_to_product_type_group[i] = new LinkedHashMap();
	
		List<String> lines  = FileUtils.readLines(new File(path));
		
		/*
		
		UNIQUE-ID - 1.14.13.78-RXN
TYPES - Chemical-Reactions
TYPES - Small-Molecule-Reactions
ATOM-MAPPINGS - (:NO-HYDROGEN-ENCODING (0 1 12 2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 20) (((ENT-KAUR-16-EN-19-OL 0 20) (NADP 21 68) (WATER 69 69)) ((CPD1F-128 0 19) (NADPH 20 67) (OXYGEN-MOLECULE 68 69))))
CREDITS - PMN
CREDITS - nilo-poyanco
ENZYMATIC-REACTION - ENZRXNN7F-45960
ENZYMATIC-REACTION - ENZRXNN7F-44282
ENZYMATIC-REACTION - ENZRXNN7F-43191
ENZYMATIC-REACTION - ENZRXNN7F-42486
ENZYMATIC-REACTION - ENZRXNN7F-41989
ENZYMATIC-REACTION - ENZRXNN7F-41925
GIBBS-0 - -85.46481    
IN-PATHWAY - PWY-5034
IN-PATHWAY - RXN-13190
INSTANCE-NAME-TEMPLATE - RXN-*
LEFT - CPD1F-128
LEFT - OXYGEN-MOLECULE
LEFT - Red-NADPH-Hemoprotein-Reductases
ORPHAN? - :NO
PHYSIOLOGICALLY-RELEVANT? - T
REACTION-BALANCE-STATUS - :BALANCED
REACTION-DIRECTION - PHYSIOL-LEFT-TO-RIGHT
RIGHT - ENT-KAUR-16-EN-19-OL
RIGHT - WATER
RIGHT - Ox-NADPH-Hemoprotein-Reductases
//SystemATIC-NAME - <em>ent</em>-kaur-16-ene,NADPH:oxygen oxidoreductase (hydroxylating)


		*/
		
		
		Boolean[] rxn_type_temp = new Boolean[rxn_type.length];
		for(int q= 0; q<rxn_type.length; q++){
		
			rxn_type_temp[q] = new Boolean(false);
		
		}
		
		String pwy_list_temp = "";
		
		
		Boolean[] formula_ele_sub_temp = new Boolean[formula_ele.length];
		for(int q= 0; q<formula_ele.length; q++){
		
			formula_ele_sub_temp[q] = new Boolean(false);
		
		}
		
		Boolean[] struc_group_sub_temp = new Boolean[struc_group.length];
		for(int q= 0; q< struc_group.length; q++){
		
			struc_group_sub_temp[q] = new Boolean(false);
		
		}
		
		Boolean[] comp_type_group_sub_temp = new Boolean[comp_type_group.length];
		for(int q= 0; q <comp_type_group.length; q++){
		
			comp_type_group_sub_temp[q] = new Boolean(false);
		
		}
		
		Boolean[] formula_ele_prod_temp = new Boolean[formula_ele.length];
		for(int q= 0; q <formula_ele.length; q++){
		
			formula_ele_prod_temp[q] = new Boolean(false);
		
		}
		
		Boolean[] struc_group_prod_temp = new Boolean[struc_group.length];
		for(int q= 0; q <struc_group.length; q++){
		
			struc_group_prod_temp[q] = new Boolean(false);
		
		}
		
		Boolean[] comp_type_group_prod_temp = new Boolean[comp_type_group.length];
		for(int q= 0; q <comp_type_group.length; q++){
		
			comp_type_group_prod_temp[q] = new Boolean(false);
		
		}
		LinkedList<String> left_id = new LinkedList();
		LinkedList<String> right_id = new LinkedList();
		//UNIQUE-ID
		//TYPES
		//ENZYMATIC-REACTION
		//IN-PATHWAY - PWY
		//LEFT -
		//RIGHT -
		
		int count = 0;
		LinkedList<String> enz_rxn_id = new LinkedList();

		for(int q= 0; q < lines.size(); q++){
		
			String one_line = lines.get(q).trim().toUpperCase();
			
			String id = "";
			
			String type = "";
			String enz_rxn = "";
			String left = "";
			String right = "";
			
			if(!one_line.equals("") && !one_line.startsWith("#")){
			
			
				if(one_line.contains("UNIQUE-ID")){
				
					if(count != 0){
					
						//System.out.println("count != 0:" + count + ":" + enz_rxn_id.size());
						//put into map with ezm_rxn id
						
						for(int c = 0; c < enz_rxn_id.size(); c++){
						
							enz_rxn = enz_rxn_id.get(c);
							if(enzrxn_id_to_rxn_type[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_rxn_type[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(rxn_type_temp[g]){
										exist[g] = rxn_type_temp[g];
									}
								}
								
								enzrxn_id_to_rxn_type[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_rxn_type[i].put(enz_rxn,rxn_type_temp);
								
							}
								enzrxn_id_to_pwy[i].put(enz_rxn,pwy_list_temp);
								
														if(enzrxn_id_to_substrate_form_ele[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_substrate_form_ele[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(formula_ele_sub_temp[g]){
										exist[g] = formula_ele_sub_temp[g];
									}
								}
								
								enzrxn_id_to_substrate_form_ele[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_substrate_form_ele[i].put(enz_rxn,formula_ele_sub_temp);
								
							}
							if(enzrxn_id_to_substrate_struc_group[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_substrate_struc_group[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(struc_group_sub_temp[g]){
										exist[g] = struc_group_sub_temp[g];
									}
								}
								
								enzrxn_id_to_substrate_struc_group[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_substrate_struc_group[i].put(enz_rxn,struc_group_sub_temp);
								
							}
							if(enzrxn_id_to_substrate_type_group[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_substrate_type_group[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(comp_type_group_sub_temp[g]){
										exist[g] = comp_type_group_sub_temp[g];
									}
								}
								
								enzrxn_id_to_substrate_type_group[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_substrate_type_group[i].put(enz_rxn,comp_type_group_sub_temp);
								
							}
							if(enzrxn_id_to_product_form_ele[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_product_form_ele[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(formula_ele_prod_temp[g]){
										exist[g] = formula_ele_prod_temp[g];
									}
								}
								
								enzrxn_id_to_product_form_ele[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_product_form_ele[i].put(enz_rxn,formula_ele_prod_temp);
								
							}
							if(enzrxn_id_to_product_struc_group[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_product_struc_group[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(struc_group_prod_temp[g]){
										exist[g] = struc_group_prod_temp[g];
									}
								}
								
								enzrxn_id_to_product_struc_group[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_product_struc_group[i].put(enz_rxn,struc_group_prod_temp);
								
							}
							if(enzrxn_id_to_product_type_group[i].containsKey(enz_rxn)){

								Boolean[] exist = enzrxn_id_to_product_type_group[i].get(enz_rxn);

								for(int g  = 0; g < exist.length; g++){

									if(comp_type_group_prod_temp[g]){
										exist[g] = comp_type_group_prod_temp[g];
									}
								}
								
								enzrxn_id_to_product_type_group[i].put(enz_rxn,exist);
							}else{
								enzrxn_id_to_product_type_group[i].put(enz_rxn,comp_type_group_prod_temp);
								
							}
							
							enzrxn_id_to_substrate_compound[i].put(enz_rxn,left_id);
							enzrxn_id_to_product_compound[i].put(enz_rxn,right_id);
						}
						
						//reset
						rxn_type_temp = new Boolean[rxn_type.length];
						for(int e= 0; e<rxn_type.length; e++){
						
							rxn_type_temp[e] = new Boolean(false);
						
						}
						
						pwy_list_temp = "";
						
						formula_ele_sub_temp = new Boolean[formula_ele.length];
						for(int e= 0; e<formula_ele.length; e++){
						
							formula_ele_sub_temp[e] = new Boolean(false);
						
						}
						
						struc_group_sub_temp = new Boolean[struc_group.length];
						for(int e= 0; e<struc_group.length; e++){
						
							struc_group_sub_temp[e] = new Boolean(false);
						
						}
						
						comp_type_group_sub_temp = new Boolean[comp_type_group.length];
						for(int e= 0; e<comp_type_group.length; e++){
						
							comp_type_group_sub_temp[e] = new Boolean(false);
						
						}
						
						formula_ele_prod_temp = new Boolean[formula_ele.length];
						for(int e= 0; e<formula_ele.length; e++){
						
							formula_ele_prod_temp[e] = new Boolean(false);
						
						}
						
						struc_group_prod_temp = new Boolean[struc_group.length];
						for(int e= 0; e<struc_group.length; e++){
						
							struc_group_prod_temp[e] = new Boolean(false);
						
						}
						
						comp_type_group_prod_temp = new Boolean[comp_type_group.length];
						for(int e= 0; e<comp_type_group.length; e++){
						
							comp_type_group_prod_temp[e] = new Boolean(false);
						
						}
						left_id = new LinkedList();
					        right_id = new LinkedList();
						enz_rxn_id = new LinkedList();
					}//count
				
					String[] space_split = space_pattern.split(one_line);
					id = space_split[space_split.length-1];
					//System.out.println("id:" + id);
					count = count +1;
					//System.out.println("count:" + count);
				
				}else if(one_line.contains("TYPES ")){
				
					String[] space_split = space_pattern.split(one_line);
					type = space_split[space_split.length-1];
					//System.out.println("rxn type:" + type);
					for(int z = 0; z < rxn_type.length; z++){
					
						if(type.equals(rxn_type[z].toUpperCase())){
						
							rxn_type_temp[z] = new Boolean(true);
							break;
						}
					}
					
				
				
				}else if(one_line.contains("ENZYMATIC-REACTION ")){
				
					String[] space_split = space_pattern.split(one_line);
					enz_rxn = space_split[space_split.length-1].trim();
					//System.out.println("enz_rxn:" + enz_rxn);
					enz_rxn_id.add(enz_rxn);
		
				
				}else if(one_line.contains("IN-PATHWAY - PWY")){
				
					String[] space_split = space_pattern.split(one_line);
					pwy_list_temp = space_split[space_split.length-1].trim();
					//System.out.println("1205 pwy:" + pwy_list_temp);
					
				
				
				}else if(one_line.contains("LEFT -")){
				
					String[] space_split = space_pattern.split(one_line);
					left = space_split[space_split.length-1].trim();
					left_id.add(left);
					//System.out.println("left:" + left);
					Boolean[] left_form = compound_id_to_formula_ele[i].get(left);

					if(left_form != null){
						for(int z = 0; z < formula_ele.length; z++){
						
							if(left_form[z]){
							
								formula_ele_sub_temp[z] = new Boolean(true);
								//System.out.println("formula_ele_sub_temp[z]:"+formula_ele_sub_temp[z]);
							
							}
							
							
						}
					}
					
					String cur_struc = compound_id_to_struc_group[i].get(left);

					if(cur_struc != null){
						for(int z = 0; z < struc_group.length; z++){
						
							if(cur_struc.contains(struc_group[z])){
								
								struc_group_sub_temp[z] = new Boolean(true);		
								//System.out.println("struc_group_sub_temp[z]:"+struc_group_sub_temp[z]);
							
							}
							
							
						}
					}
					
					String cur_type_group = compound_id_to_type_group[i].get(left);

					if(cur_type_group != null){
						for(int z = 0; z < comp_type_group.length; z++){
						
							if(cur_type_group.contains(comp_type_group[z])){
								
								comp_type_group_sub_temp[z] = new Boolean(true);
								//System.out.println("comp_type_group_sub_temp[z] :"+comp_type_group_sub_temp[z] );
							}
							
							
						}
					}
				
				}else if(one_line.contains("RIGHT -")){
				
					String[] space_split = space_pattern.split(one_line);
					right = space_split[space_split.length-1];
					right_id.add(right);
					
					Boolean[] right_form = compound_id_to_formula_ele[i].get(right);

					if(right_form != null){
						for(int z = 0; z < formula_ele.length; z++){
						
							if(right_form[z]){
							
								
							
								formula_ele_prod_temp[z] = new Boolean(true);
								//System.out.println("formula_ele_prod_temp[z]:"+formula_ele_prod_temp[z] );
							}
							
							
						}
					}
					
					String cur_right_struc = compound_id_to_struc_group[i].get(right);
					
					if(cur_right_struc != null){
						for(int z = 0; z < struc_group.length; z++){
						
							if(cur_right_struc.contains(struc_group[z])){
								
								struc_group_prod_temp[z] = new Boolean(true);
								//System.out.println("struc_group_prod_temp[z]:"+struc_group_prod_temp[z] );
								
							}
							
							
						}
					}
					
					String cur_right_type_group = compound_id_to_type_group[i].get(right);

					if(cur_right_type_group != null){
						for(int z = 0; z < comp_type_group.length; z++){
						
							if(cur_right_type_group.contains(comp_type_group[z])){
							
								comp_type_group_prod_temp[z] = new Boolean(true);
								//System.out.println("comp_type_group_prod_temp[z]:"+comp_type_group_prod_temp[z]);
							
							}
							
							
						}
					}
				
				}

				
			}
		
		}//q
		
		
		
		LinkedList<String> keys = new LinkedList<String>(enzrxn_id_to_rxn_type[i].keySet());
		////System.out.println(i + ":" + "keys.size():"+ keys.size());
		for(int j = 0; j < keys.size(); j++){
		
			String rxn_id = keys.get(j);
			Boolean[] rxn_type = enzrxn_id_to_rxn_type[0].get(rxn_id);
			String pwy  = enzrxn_id_to_pwy[0].get(rxn_id);
			LinkedList<String> sub_com = enzrxn_id_to_substrate_compound[0].get(rxn_id);
			Boolean[] sub_form = enzrxn_id_to_substrate_form_ele[0].get(rxn_id);
			Boolean[] sub_str = enzrxn_id_to_substrate_struc_group[0].get(rxn_id);
			Boolean[] sub_type = enzrxn_id_to_substrate_type_group[0].get(rxn_id);
			LinkedList<String> pro_com = enzrxn_id_to_product_compound[0].get(rxn_id);
			Boolean[] pro_form = enzrxn_id_to_product_form_ele[0].get(rxn_id);
			Boolean[] pro_str = enzrxn_id_to_product_struc_group[0].get(rxn_id);
			Boolean[] pro_type = enzrxn_id_to_product_type_group[0].get(rxn_id);
			
			//System.out.println("error rxn_id:"+rxn_id + ":" + rxn_id);
			//System.out.println("error rxn_type:"+rxn_id + ":" +Arrays.toString(rxn_type));
			//System.out.println("error pwy:"+rxn_id + ":" +pwy);
			//System.out.println("error sub_form:"+rxn_id + ":" + Arrays.toString(sub_form));
			//System.out.println("error sub_str:"+rxn_id + ":" + Arrays.toString(sub_str));
			//System.out.println("error sub_type:"+rxn_id + ":" + Arrays.toString(sub_type));
			//System.out.println("error pro_form:"+rxn_id + ":" + Arrays.toString(pro_form));
			
			//System.out.println("error pro_str:"+rxn_id + ":" + Arrays.toString(pro_str));
			//System.out.println("error pro_type:"+rxn_id + ":" +Arrays.toString(pro_type));
			//System.out.println("error sub_com:" +rxn_id + ":" + Arrays.toString(sub_com.toArray()));
			//System.out.println("error pro_com:" +rxn_id + ":" + Arrays.toString(pro_com.toArray()));
			
		}

	}

}


//genes.dat
static LinkedHashMap<String,String>[] gene_name_to_pmn_acc;
static LinkedHashMap<String,String>[] product_to_pmn_acc;
static LinkedHashMap<String, LinkedList<String>>[] pmn_acc_to_product;

void  map_gene_name_to_pmn_acc()  throws IOException{

	//System.out.println("map_gene_name_to_pmn_acc():" );

	gene_name_to_pmn_acc = new LinkedHashMap[species.length];
	product_to_pmn_acc = new LinkedHashMap[species.length];
	pmn_acc_to_product = new LinkedHashMap[species.length];
	
	String[] species_file_path = {"big/pmn/aracyc/17.1.2/data/"};
	//,"big/pmn/brachypodiumcyc/8.0.0/data/","big/pmn/watermeloncyc/1.0.2/data/","big/pmn/csativa_pkcyc/3.0.2/data/","big/pmn/graimondiicyc/3.0.2/data/","big/pmn/mdomesticacyc/3.0.2/data/","big/pmn/oryzacyc/7.2.1/data/","big/pmn/commonbeancyc/3.0.2/data/","big/pmn/poplarcyc/13.0.0/data/","big/pmn/ppersicacyc/3.0.2/data/","big/pmn/setariacyc/7.1.1/data/","big/pmn/tomatocyc/6.0.0/data/","big/pmn/potatocyc/7.0.0/data/","big/pmn/sorghumbicolorcyc/7.1.1/data/","big/pmn/grapecyc/10.0.0/data/","big/pmn/corncyc/11.0.0/data/"};

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	for(int i = 0; i < species_file_path.length; i++){

		String path = species_file_path[i] + "genes.dat";
		gene_name_to_pmn_acc[i] = new LinkedHashMap();
		product_to_pmn_acc[i] = new LinkedHashMap();	
		pmn_acc_to_product[i] = new LinkedHashMap();	
	
		List<String> lines  = FileUtils.readLines(new File(path));
		
		String id = "";
		String acc = "";
		String prot_prod = "";


		/*
		UNIQUE-ID - OS10G0195600
TYPES - ORFs
COMMON-NAME - Os10g0195600
ACCESSION-1 - Os10g0195600
INSTANCE-NAME-TEMPLATE - G-*
PRODUCT - OS10G0195600-MONOMER
UNMAPPED-COMPONENT-OF - ORYZA-SATIVA
		*/
					
		for(int k = 0; k < lines.size(); k++){
					
			String one_line = lines.get(k).toUpperCase().trim();
			if(!one_line.equals("") && !one_line.startsWith("#")){
			
				
				if(one_line.contains("UNIQUE-ID ")){
				
					String[] space_split = space_pattern.split(one_line.trim());
					id = space_split[space_split.length-1].trim();
					//System.out.println("1416:" + id);
				}else if(one_line.contains("ACCESSION-1 -")){
				
					String[] space_split = space_pattern.split(one_line);
					acc = space_split[space_split.length-1].trim();
					//System.out.println("1416:" + acc);
					gene_name_to_pmn_acc[i].put(id,acc);
				}else if(one_line.contains("PRODUCT - ")){
				
					String[] space_split = space_pattern.split(one_line);
					prot_prod = space_split[space_split.length-1].trim();
					//System.out.println("1416:" + prot_prod);
					product_to_pmn_acc[i].put(prot_prod,acc);
					
					if(pmn_acc_to_product[i].containsKey(acc)){
						LinkedList<String> exist = pmn_acc_to_product[i].get(acc);
						exist.add(prot_prod);
						pmn_acc_to_product[i].put(acc,exist);
					}else{
						LinkedList<String> exist = new LinkedList();
						exist.add(prot_prod);
						
						pmn_acc_to_product[i].put(acc,exist);
					
					}
				}
			}
					
			
		}
		
	}
}

//enzrxns.dat

//UNIQUE-ID
//ENZYME - 

static LinkedHashMap<String,String>[] product_to_enzrxn_id;
static LinkedHashMap<String,LinkedList<String>>[] enzrxn_id_to_product;

void  map_product_to_enzrxn_id()  throws IOException{

	//System.out.println("map_product_to_enzrxn_id()" );

	product_to_enzrxn_id = new LinkedHashMap[species.length];
	enzrxn_id_to_product = new LinkedHashMap[species.length];
	
	String[] species_file_path = {"big/pmn/aracyc/17.1.2/data/"};
	
	//,"big/pmn/brachypodiumcyc/8.0.0/data/","big/pmn/watermeloncyc/1.0.2/data/","big/pmn/csativa_pkcyc/3.0.2/data/","big/pmn/graimondiicyc/3.0.2/data/","big/pmn/mdomesticacyc/3.0.2/data/","big/pmn/oryzacyc/7.2.1/data/","big/pmn/commonbeancyc/3.0.2/data/","big/pmn/poplarcyc/13.0.0/data/","big/pmn/ppersicacyc/3.0.2/data/","big/pmn/setariacyc/7.1.1/data/","big/pmn/tomatocyc/6.0.0/data/","big/pmn/potatocyc/7.0.0/data/","big/pmn/sorghumbicolorcyc/7.1.1/data/","big/pmn/grapecyc/10.0.0/data/","big/pmn/corncyc/11.0.0/data/"};

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	for(int i = 0; i < species_file_path.length; i++){

		String path = species_file_path[i] + "enzrxns.dat";
		product_to_enzrxn_id[i] = new LinkedHashMap();
		enzrxn_id_to_product[i] = new LinkedHashMap();
		List<String> lines  = FileUtils.readLines(new File(path));
		
		String id = "";
		String enz = "";
		

		/*
		UNIQUE-ID - ENZRXNN7F-50378
TYPES - Enzymatic-Reactions
COMMON-NAME - acyl-acp:<i>sn</i>-glycerol-3-phosphate 1-<i>O</i>-acyltransferase
BASIS-FOR-ASSIGNMENT - :MANUAL
CITATIONS - E2P2PMN2019:EV-COMP-AINF
CREDITS - pmngroup
ENZYME - OS11G0679700-MONOMER
INSTANCE-NAME-TEMPLATE - ENZRXN-*
PHYSIOLOGICALLY-RELEVANT? - T
REACTION - RXN-10462
		*/
					
		for(int k = 0; k < lines.size(); k++){
					
			String one_line = lines.get(k).toUpperCase();
			if(!one_line.equals("") && !one_line.startsWith("#")){
			
			
				if(one_line.contains("UNIQUE-ID ")){
				
					String[] space_split = space_pattern.split(one_line);
					id = space_split[space_split.length-1];
					//System.out.println("enzrxn id:" + id);
				}else if(one_line.contains("ENZYME - ")){
				
					String[] space_split = space_pattern.split(one_line);
					enz = space_split[space_split.length-1];
					//enz = enz.substring(0, enz.indexOf("-"));
					//System.out.println("product enz:" + enz);
					//String pmn_acc = product_to_pmn_acc[i].get(enz);
					////System.out.println("pmn_acc:" +pmn_acc);
					//if(pmn_acc != null){
						////System.out.println("1511:" +pmn_acc + ":" + id);
						product_to_enzrxn_id[i].put(enz,id);
						if(enzrxn_id_to_product[i].containsKey(id)){

								LinkedList<String> exist = enzrxn_id_to_product[i].get(id);

								if(!exist.contains(enz)){
									exist.add(enz);
								}
		
								
								enzrxn_id_to_product[i].put(id,exist);
							}else{

								LinkedList<String> exist = new LinkedList();
								exist.add(enz);
								enzrxn_id_to_product[i].put(id,exist);
								
							}

					
					//}
					/*if(!enz.matches("^OS\\d+G\\d+.*")){

						//enz = gene_name_to_pmn_acc.get(enz);
						if(enz == null){

							enz = product_to_pmn_acc.get(enz);

							if(enz != null){
								if(enz.matches("^OS\\d+G\\d+.*")){

									product_to_enzrxn_id[i].put(enz,id);

								}

							}
						}else{

							if(enz.matches("^OS\\d+G\\d+.*")){

								product_to_enzrxn_id[i].put(enz,id);

							}

							

						}

					}else{

					
						product_to_enzrxn_id[i].put(enz,id);
					}
					*/
				}
			}
					
			
		}
		
	}
}//method

static LinkedHashMap<String,LinkedList<String>>[] pwy_to_acc;
static LinkedList<String>[] pwy_list_to_acc;

void  from_pwy_to_acc()  throws IOException{ 

	pwy_list_to_acc = new LinkedList[pwy_list.length];
	pwy_to_acc = new LinkedHashMap[species.length];
	for(int i = 0; i < pwy_list.length; i++){

		pwy_list_to_acc[i] = new LinkedList();
		
	}
	
	for(int i = 0; i < species.length; i++){
	
		pwy_to_acc[i] = new LinkedHashMap();
	}

	for(int i = 0; i < species.length; i++){

		LinkedList<String> keys = new LinkedList<String>(enzrxn_id_to_pwy[i].keySet());
		////System.out.println(i + ":" + "keys.size():"+ keys.size());
		for(int j = 0; j < keys.size(); j++){

			String key = keys.get(j);
			////System.out.println(i + ":" + "key:"+ key);
			String pwy_array = enzrxn_id_to_pwy[i].get(key);
			for(int k = 0; k < pwy_list.length; k++){

				if(pwy_list[k].equals(pwy_array)){

					LinkedList<String> acc_list = enzrxn_id_to_product[i].get(key);
					if(acc_list != null){
						////System.out.println("acc_list:"+ Arrays.toString(acc_list.toArray()));
						pwy_list_to_acc[k].addAll(acc_list);
					}
					
				}
			}
		}//j

		for(int j = 0; j < pwy_list.length; j++){

			String cur_pwy = pwy_list[j];
			LinkedList<String> val = pwy_list_to_acc[j];
			Set<String> pfams_h_inv = new LinkedHashSet<>();
			pfams_h_inv.addAll((List)val);
										
			// Clear the list
			((List)val).clear();
																  
			// add the elements of set
			// with no duplicates to the list
			((List)val).addAll(pfams_h_inv);
			pwy_to_acc[i].put(cur_pwy,val);
			//System.out.println("pwy_to_gene:"+ cur_pwy + ":" + Arrays.toString(val.toArray()));
			
		}

			
	}//i
					
}//method


static String[] super_path = {
"Acetyl-CoA-Biosynthesis","Allantoin-degradation","Amino-Acid-Biosynthesis","Ammonia-Assimilation","ANTHOCYANIN-SYN","ARGININE-SYN","ASPARAGINE-SYN","Benzoate-Biosynthesis","C40-Carotenoids-Biosynthesis","Chorismate-Biosynthesis","Citrulline-Biosynthesis","CoA-Biosynthesis","COUMARIN-SYN","D-Amino-Acid-Degradation","Energy-Metabolism","Fatty-acid-biosynthesis","FATTY-ACID-DERIVATIVE-SYN","Folate-Biosynthesis","GGPP-Biosynthesis","GIBBERELLIN-SYN","Guanosine-Nucleotides-Degradation","Heme-b-Biosynthesis","HEME-SYN","Indole-3-Acetate-Inactivation","Interconversion","ISOLEUCINE-SYN","Lipid-Biosynthesis","Methionine-Salvage","Pentose-Phosphate-Cycle","PhosphatidylglycerolBiosynthesis","Phospholipid-Biosynthesis","Photosynthesis","Phylloquinone-Biosynthesis","Plastoquinone-Biosynthesis","Polyamine-Biosynthesis","Purine-Degradation","Purine-Nucleotide-De-Novo-Biosynthesis","Purine-Nucleotides-Salvage","Putrescine-Biosynthesis","Pyrimid-Deoxyribonucleot-De-Novo-Biosyn","Pyrimidine-Nucleotide-Salvage","Pyrimid-Ribonucleot-De-Novo-Biosyn","REACTIVE-OXYGEN-SPECIES-DEGRADATION","Reductants","Seleno-Amino-Acid-Detoxification","Sterol-Biosynthesis","Sucrose-Biosynthesis","SUCROSE-DEG"

};


static LinkedHashMap<String,String>[] rxn_id_to_enzrxn_id;
void map_rxn_id_to_enzrxn_id() throws IOException{

	//System.out.println("map_rxn_id_to_enzrxn_id:" );
	rxn_id_to_enzrxn_id = new LinkedHashMap[species.length];
	rxn_id_to_enzrxn_id[0] = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	

	String path = "big/pmn/aracyc/17.1.2/data/reactions.dat";
	List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){

  		String one_line = lines.get(q).trim().toUpperCase();
  		if(one_line.contains("UNIQUE-ID - ")){
  		
  			String rxn_id = one_line.replaceFirst("UNIQUE-ID - ","");
  			
  			for(int w= q; w < lines.size(); w++){
  			
  				String next_line = lines.get(w).trim().toUpperCase();
  				
  				if(next_line.contains("UNIQUE-ID - ")){
  					break;
  				}
  			
  				if(next_line.contains("ENZYMATIC-REACTION - ")){
  				
  					String enzrxn_id = next_line.replaceFirst("ENZYMATIC-REACTION - ","");
  					
  					System.out.println(rxn_id + ":" + enzrxn_id);
  					rxn_id_to_enzrxn_id[0].put(rxn_id,enzrxn_id);
  				}
  			}
  		}
		
	}

}

static LinkedHashMap<String,LinkedList<String>>[] map_pathway_to_super_pathway;

void retrieve_super_pathway() throws IOException{

	//System.out.println("retrieve_plantcyc_comp_type:" );
	map_pathway_to_super_pathway = new LinkedHashMap[species.length];
	
	for(int i = 0; i < species.length; i++){
	
		map_pathway_to_super_pathway[i] = new LinkedHashMap();
	}

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern dash_pattern = Pattern.compile("-");

	String path = "big/pmn/aracyc/17.1.2/data/pathways.dat";

	
	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		if(one_line.contains("UNIQUE-ID - ")){

			//System.out.println("498 one_line:" + one_line);
			String id = one_line.replaceFirst("UNIQUE-ID - ","");

			for(int r = q; r < lines.size(); r++){
			
				String next_line = lines.get(r).trim();
			
				if(next_line.contains("UNIQUE-ID - ")){
				
					break;
				}

				
				if(next_line.contains("SUPER-PATHWAYS - ")){
				
					String super_id = next_line.replaceFirst("SUPER-PATHWAYS - ","");
				
					if(map_pathway_to_super_pathway[0].containsKey(id)){

								LinkedList<String> exist = map_pathway_to_super_pathway[0].get(id);

								if(!exist.contains(super_id)){
								
									exist.add(super_id);
								
								}
								map_pathway_to_super_pathway[0].put(id,exist);
								
							}else{
							
								LinkedList<String> exist = new LinkedList();
								exist.add(super_id);
								map_pathway_to_super_pathway[0].put(id,exist);
								
							}
				}
			}
		}
	}

				
}//method




}//inner class pathway analysis

static class multiomics_asso_module{

static String[] type = {"Promoter","Enhancer"};

static String[] study_16_factor = {"arsenic","cadmium","cobalt","copper","potassium","magnesium","manganese","molybdenum","Nickel","phosphorus","rubidium","selenium","sulfur","zinc"};

static String[] study_16_pheno ={"As75","Cd111","Co59","Cu65","K39","Mg25","Mn55","Mo98","Ni60","P31","Rb85","S34","Se82","Zn66"};

static String[] epi = {"CG","CHG","CHH"};

static String[] main_path = {"ALACAT2-PWY","ALANINE-DEG3-PWY","ALANINE-SYN2-PWY","ARGASEDEG-PWY","ARGDEG-V-PWY","ARG-PRO-PWY","ARGSPECAT-PWY","ARGSYNBSUB-PWY","ARGSYN-PWY","ASPARTATE-DEG1-PWY","ASPARTATESYN-PWY","ASPASN-ARA-PWY","ASPSYNII-PWY","BSUBPOLYAMSYN-PWY","CALVIN-PWY","CHLOROPHYLL-SYN","CITRULBIO-PWY","CITRULLINE-DEG-PWY","CYANCAT-PWY","CYSTSYN-PWY","ETHYL-PWY","FASYN-ELONG-PWY","GLNSYN-PWY","GLUGLNSYN-PWY","GLUTAMINDEG-PWY","GLUTATHIONESYN-PWY","GLUTORN-PWY","GLUT-REDOX-PWY","GLUTSYNIII-PWY","GLYOXYLATE-BYPASS","GLYSYN2-PWY","GLYSYN-ALA-PWY","GLYSYN-PWY","HEME-BIOSYNTHESIS-II","HOMOSERSYN-PWY","HOMOSER-THRESYN-PWY","ILEUDEG-PWY","ILEUSYN-PWY","KDOSYN-PWY","LEU-DEG2-PWY","LEUSYN-PWY","LIPAS-PWY","LYSINE-DEG2-PWY","MALATE-ASPARTATE-SHUTTLE-PWY","METHIONINE-DEG1-PWY","NAGLIPASYN-PWY","NONMEVIPP-PWY","NONOXIPENT-PWY","OXIDATIVEPENT-PWY","P401-PWY","PANTO-PWY","PLPSAL-PWY","POLYAMINSYN3-PWY","PROSYN-PWY","PROUT-PWY","PWY0-1182","PWY0-1264","PWY0-1313","PWY0-1507","PWY0-1535","PWY0-166","PWY0-662","PWY0A-6303","PWY-1001","PWY-102","PWY-1061","PWY-1081","PWY-1121","PWY-1186","PWY-1187","PWY-1422","PWY-1581","PWY-1782","PWY-181","PWY-1822","PWY-1881","PWY1F-353","PWY1F-467","PWY1F-823","PWY-2","PWY-2021","PWY-2161","PWY-2301","PWY-2541","PWY-2582","PWY-2602","PWY-2681","PWY-2781","PWY-282","PWY-283","PWY-2841","PWY-2881","PWY-2901","PWY-2902","PWY-3001","PWY-3041","PWY-3301","PWY-3341","PWY-3385","PWY-3461","PWY-3462","PWY-361","PWY-3781","PWY-3821","PWY-3841","PWY-3982","PWY-4","PWY-40","PWY-401","PWY-4041","PWY-4101","PWY-4203","PWY-4261","PWY-43","PWY-4302","PWY-4341","PWY-4361","PWY-4381","PWY-4661","PWY-4821","PWY-4861","PWY490-4","PWY-4984","PWY4FS-11","PWY4FS-12","PWY4FS-13","PWY4FS-2","PWY4FS-3","PWY4FS-4","PWY4FS-6","PWY4FS-7","PWY4FS-8","PWY-5027","PWY-5032","PWY-5034","PWY-5035","PWY-5041","PWY-5046","PWY-5049","PWY-5060","PWY-5063","PWY-5064","PWY-5068","PWY-5070","PWY-5083","PWY-5084","PWY-5098","PWY-5107","PWY-5120","PWY-5121","PWY-5123","PWY-5125","PWY-5129","PWY-5136","PWY-5137","PWY-5138","PWY-5143","PWY-5147","PWY-5152","PWY-5172","PWY-5173","PWY-5175","PWY-5194","PWY-5203","PWY-5267","PWY-5268","PWY-5269","PWY-5271","PWY-5285","PWY-5287","PWY-5313","PWY-5326","PWY-5337","PWY-5340","PWY-5377","PWY-5391","PWY-5410","PWY-5434","PWY-5441","PWY-5453","PWY-5481","PWY-5530","PWY-561","PWY-5640","PWY-5667","PWY-5669","PWY-5670","PWY-5686","PWY-5687","PWY-5690","PWY-5691","PWY-5692","PWY-5697","PWY-5698","PWY-5704","PWY-5723","PWY-5725","PWY-5788","PWY-5791","PWY-5800","PWY-5807","PWY-581","PWY-5837","PWY-5859","PWY-5863","PWY-5871","PWY-5876","PWY-5886","PWY-5905","PWY-5910","PWY-5921","PWY-5934","PWY-5936","PWY-5943","PWY-5945","PWY-5946","PWY-5947","PWY-5971","PWY-5973","PWY-5980","PWY-5986","PWY-5989","PWY-5992","PWY-5997","PWY-5998","PWY-6","PWY-6005","PWY-6007","PWY-6008","PWY-6012","PWY-6019","PWY-6039","PWY-6064","PWY-6066","PWY-6121","PWY-6122","PWY-6124","PWY-6132","PWY-6137","PWY-6147","PWY-6163","PWY-6164","PWY-6196","PWY-6199","PWY-621","PWY-6219","PWY-6232","PWY-6266","PWY-6287","PWY-6303","PWY-6305","PWY-6364","PWY-641","PWY-6424","PWY-6441","PWY-6442","PWY-6443","PWY-6457","PWY-6473","PWY-6477","PWY-6498","PWY-6502","PWY-6527","PWY-6535","PWY-6544","PWY-6546","PWY-6549","PWY-6556","PWY-6596","PWY-6599","PWY-66","PWY-6605","PWY-6606","PWY-6607","PWY-6613","PWY-6619","PWY66-423","PWY-6663","PWY-6668","PWY-6673","PWY-6707","PWY-6730","PWY-6733","PWY-6736","PWY-6754","PWY-6756","PWY-6773","PWY-6786","PWY-6792","PWY-6799","PWY-6802","PWY-6804","PWY-6806","PWY-6825","PWY-6837","PWY-6909","PWY-6922","PWY-6930","PWY-6931","PWY-6932","PWY-6936","PWY-695","PWY-6952","PWY-6959","PWY-6963","PWY-6964","PWY-699","PWY-7036","PWY-7047","PWY-7048","PWY-7060","PWY-7061","PWY-7067","PWY-7068","PWY-7071","PWY-7101","PWY-7140","PWY-7170","PWY-7176","PWY-7177","PWY-7183","PWY-7184","PWY-7185","PWY-7187","PWY-7196","PWY-7197","PWY-7205","PWY-7208","PWY-7219","PWY-7221","PWY-7224","PWY-7226","PWY-7227","PWY-7267","PWY-7270","PWY-7325","PWY-7343","PWY-7344","PWY-735","PWY-7388","PWY-7398","PWY-7416","PWY-7417","PWY-7436","PWY-7450","PWY-7528","PWY-7560","PWY-7590","PWY-7618","PWY-7640","PWY-782","PWY-7856","PWY-7859","PWY-7909","PWY-801","PWY-83","PWY-861","PWY-882","PWY-922","PWYDQC-4","PWYQT-4429","PWYQT-4432","PWYQT-4433","PWYQT-4445","PWYQT-4450","PWYQT-4470","PWYQT-4471","PWYQT-4472","PWYQT-4473","PWYQT-4474","PWYQT-4475","PWYQT-4476","PWYQT-4477","PWYQT-4481","PWYQT-62","PYRIDNUCSYN-PWY","PYRUVDEHYD-PWY","RIBOSYN2-PWY","SAM-PWY","SULFMETII-PWY","THIOREDOX-PWY","THISYNARA-PWY","THRESYN-PWY","TRESYN-PWY","TYRFUMCAT-PWY","VALDEG-PWY","VALSYN-PWY","XYLCAT-PWY"};

static LinkedHashMap<String,Boolean[]> gene_to_pheno;
static LinkedHashMap<String,LinkedList<Double>[][][]> gene_to_all_ewas;
static LinkedHashMap<String,Double[][][]> gene_to_avr_ewas;
	
void map_pheno_to_path_study16() throws IOException{

	System.out.println("map_pheno_to_path_main_deg():" );

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern digit_pattern = Pattern.compile("\\d");

	gene_to_pheno  = new LinkedHashMap();
	gene_to_all_ewas  = new LinkedHashMap();
	gene_to_avr_ewas  = new LinkedHashMap();

	String path = "data/out_rewrite_ewas_enhancer_promoter";

	//
	//1	77382	77383	K39	CHH	-0.380934726820547	enhancer	AT1G01140

	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
						
		String line = lines.get(q).trim();
		String[] tab_split = tab_pattern.split(line);
		String gene = tab_split[7].trim();
		String cur_pheno = tab_split[3].trim();
		String cur_epi = tab_split[4].trim();
		Double cur_val = new Double(tab_split[5].trim());
		String cur_type = tab_split[6].trim();

		int type_int = -1;
		if(cur_type.equals("Promoter")){

			type_int = 0;
		}else{

			type_int = 1;
		}

		int epi_int = -1;

		for(int i = 0; i < epi.length; i++){
	
			if(cur_epi.equals(epi[i])){

				epi_int = i;
				break;
			}
		}

		int pheno_int = -1;

		for(int i = 0; i < study_16_pheno.length; i++){
	
			if(cur_pheno.equals(study_16_pheno[i])){

				pheno_int = i;
				break;
			}
		}

		if(gene_to_pheno.containsKey(gene)){
			
			Boolean[] exist = gene_to_pheno.get(gene);
			for(int r = 0; r < study_16_pheno.length; r++){
			
				if(cur_pheno.equals(study_16_pheno[r])){
					exist[r] = new Boolean(true);
				}
			}
			gene_to_pheno.put(gene,exist);
				
		}else{
			Boolean[] exist = new Boolean[study_16_pheno.length];
			for(int r = 0; r < study_16_pheno.length; r++){
			
				if(cur_pheno.equals(study_16_pheno[r])){
					exist[r] = new Boolean(true);
				}else{

					exist[r] = new Boolean(false);
				}
			}
			gene_to_pheno.put(gene,exist);
		}

		if(gene_to_all_ewas.containsKey(gene)){
			
			LinkedList<Double>[][][] exist = gene_to_all_ewas.get(gene);
			exist[type_int][pheno_int][epi_int].add(cur_val);
			
			gene_to_all_ewas.put(gene,exist);
				
		}else{
			LinkedList<Double>[][][] exist = new LinkedList[2][study_16_pheno.length][epi.length];
			for(int z = 0; z < type.length; z++){

				for(int x = 0; x < study_16_pheno.length; x++){

					for(int y = 0; y < epi.length; y++){

						exist[z][x][y] = new LinkedList();

						if(z == type_int && x == pheno_int && y == epi_int){

							exist[z][x][y].add(cur_val);
						}

					}
				}
			}
			gene_to_all_ewas.put(gene,exist);
			
		}
	}//lines

	List<String> keys_list = new LinkedList<String>(gene_to_all_ewas.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		String gene = keys_list.get(i);
		Double[][][] input = new Double[type.length][study_16_pheno.length][epi.length];

		LinkedList<Double>[][][] exist = gene_to_all_ewas.get(gene);

		for(int z = 0; z < exist.length; z++){
			for(int x = 0; x < exist[z].length; x++){
				for(int y = 0; y < exist[z][x].length; y++){

					double[] values = ArrayUtils.toPrimitive(exist[z][x][y].toArray(new Double[0]));

					if(values.length == 1){

						input[z][x][y] = values[0];
					}else{


						DescriptiveStatistics desc = new DescriptiveStatistics(values);
						input[z][x][y]=new Double(desc.getMean());
					}
				}
       			}
		}
			
			
			
		gene_to_avr_ewas.put(gene,input);
				
	}
}

//has isoform

static ArrayList<String> isoform_flag;
static LinkedHashMap<String,Integer> gene_to_isoform;
	
void get_isoform_flag() throws IOException{

	System.out.println("get_isoform_flag() :" );

	Pattern space_pattern = Pattern.compile("\\s+");
	
	isoform_flag  = new ArrayList();
	gene_to_isoform  = new LinkedHashMap();
	
	String path = "data/Araport11_GFF3_genes_transposons.current.gff.mRNA.sort.uniqc";

	/*
	      1 AT1G01010
      	6 AT1G01020
      */
	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
						
		String line = lines.get(q).trim();
		String[] tab_split = space_pattern.split(line);
		String gene = tab_split[1].trim();
		String num = tab_split[0].trim();
		
		int num_int = new Integer(num).intValue();
		
		if(num_int>=2){
		
			isoform_flag.add(gene);
			gene_to_isoform.put(gene, new Integer(num_int));
		}
		
	}
}//method

static LinkedHashMap<String,Boolean[]> gene_to_pheno_twas_pos;
static LinkedHashMap<String,Boolean[]> gene_to_pheno_twas_neg;
//study_16_pheno 

	
void map_gene_to_pheno_twas() throws IOException{

	System.out.println("map_gene_to_pheno_twas():" );

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern digit_pattern = Pattern.compile("\\d");

	gene_to_pheno_twas_pos  = new LinkedHashMap();
	gene_to_pheno_twas_neg  = new LinkedHashMap();
	
	String path = "data/twas.pheno.csv";

	//Co59,AT1G09220,-0.458974767,9.85E-34,5.49E-33,study16,-0.562125698,4.21E-53,"5.25e-52",cobalt concentration
	
	ArrayList<String> pheno_list = new ArrayList();
	
	for(int i = 0; i < study_16_pheno.length; i++){
		pheno_list.add(study_16_pheno[i]);
	
	}

	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
						
		String line = lines.get(q).trim();
		String[] tab_split = comma_pattern.split(line);
		
		if(pheno_list.contains(tab_split[0].trim())){
		
			String gene = tab_split[1].trim();
			String pearson1 = tab_split[2].trim();
			String pearson2 = tab_split[6].trim();

			
			if(new Double(pearson1) >= 0.5 || new Double(pearson2) >= 0.5){
			
				int study_type_int = -1;
				for(int i = 0; i < study_16_pheno.length; i++){
				
					if(tab_split[0].trim().equals(study_16_pheno[i])){
						study_type_int = i;
						break;
					}
				
				}
			
				if(gene_to_pheno_twas_pos.containsKey(gene)){
				
					Boolean[] exist = gene_to_pheno_twas_pos.get(gene);
					exist[study_type_int] = new Boolean(true);
					gene_to_pheno_twas_pos.put(gene,exist);
						
				}else{
					Boolean[] exist = new Boolean[study_16_pheno.length];
					for(int r = 0; r < study_16_pheno.length; r++){
					
						if(r==study_type_int){
							exist[r] = new Boolean(true);
						}else{

							exist[r] = new Boolean(false);
						}
					}
					gene_to_pheno_twas_pos.put(gene,exist);
				}
			}//if greater or less
			
			if(new Double(pearson1) <= -0.5 || new Double(pearson2) <= -0.5){
			
				int study_type_int = -1;
				for(int i = 0; i < study_16_pheno.length; i++){
				
					if(tab_split[0].trim().equals(study_16_pheno[i])){
						study_type_int = i;
						break;
					}
				
				}
			
				if(gene_to_pheno_twas_neg.containsKey(gene)){
				
					Boolean[] exist = gene_to_pheno_twas_neg.get(gene);
					exist[study_type_int] = new Boolean(true);
					gene_to_pheno_twas_neg.put(gene,exist);
						
				}else{
					Boolean[] exist = new Boolean[study_16_pheno.length];
					for(int r = 0; r < study_16_pheno.length; r++){
					
						if(r==study_type_int){
							exist[r] = new Boolean(true);
						}else{

							exist[r] = new Boolean(false);
						}
					}
					gene_to_pheno_twas_neg.put(gene,exist);
				}
			}//if greater or less
		}
	}
	
}//method

//gene to pheno twas


static LinkedHashMap<String,Boolean[]> pheno_to_pwy;
static LinkedHashMap<String,LinkedList<Double>[]> pheno_to_all_pwy;
static LinkedHashMap<String,Double[]> pheno_to_avr_pwy;

void map_pheno_to_path() throws IOException{

	System.out.println("map_pheno_to_path():" );

	pheno_to_pwy = new LinkedHashMap();
	pheno_to_all_pwy = new LinkedHashMap();
	pheno_to_avr_pwy = new LinkedHashMap();
/*
As75,PWY1F-467,0.345727717,6.25E-19,6.10E-17,study16,0.311841692,1.62E-15,"8.63e-14
",arsenic concentration

*/

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern digit_pattern = Pattern.compile("\\d");

	String path = "data/phenotypepathway.arapheno.study16.intersect";

	

	List<String> lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < lines.size(); q++){
						
		String line = lines.get(q).trim();
		line = line.replaceAll("\"","");
		String[] comma_split = comma_pattern.split(line);
		String pheno = comma_split[0].trim();
		String pwy = comma_split[1].trim();

		if(comma_split.length >=9){
			if(pheno_to_pwy.containsKey(pheno)){
				
				Boolean[] exist = pheno_to_pwy.get(pheno);
				for(int r = 0; r < main_path.length; r++){
				
					if(pwy.equals(main_path[r])){
						exist[r] = new Boolean(true);
					}
				}
				pheno_to_pwy.put(pheno,exist);
					
			}else{
				Boolean[] exist = new Boolean[main_path.length];
				for(int r = 0; r < main_path.length; r++){
				
					if(pwy.equals(main_path[r])){
						exist[r] = new Boolean(true);
					}else{

						exist[r] = new Boolean(false);
					}
				}
				pheno_to_pwy.put(pheno,exist);
			}

			if(pheno_to_all_pwy.containsKey(pheno)){
				
				LinkedList<Double>[] exist = pheno_to_all_pwy.get(pheno);
				exist[0].add(new Double(comma_split[2].trim()));
				exist[1].add(new Double(comma_split[3].trim()));
				exist[2].add(new Double(comma_split[4].trim()));
				exist[3].add(new Double(comma_split[6].trim()));
				exist[4].add(new Double(comma_split[7].trim()));
				exist[5].add(new Double(comma_split[8].trim()));
				pheno_to_all_pwy.put(pheno,exist);
					
			}else{
				LinkedList<Double>[] exist = new LinkedList[6];
				exist[0] = new LinkedList(); exist[0].add(new Double(comma_split[2].trim()));
				exist[1] = new LinkedList(); exist[1].add(new Double(comma_split[3].trim()));
				exist[2] = new LinkedList(); exist[2].add(new Double(comma_split[4].trim()));
				exist[3] = new LinkedList(); exist[3].add(new Double(comma_split[6].trim()));
				exist[4] = new LinkedList(); exist[4].add(new Double(comma_split[7].trim()));
				exist[5] = new LinkedList(); exist[5].add(new Double(comma_split[8].trim()));
				pheno_to_all_pwy.put(pheno,exist);
				
			}
		}
			

	}//lines

	List<String> keys_list = new LinkedList<String>(pheno_to_all_pwy.keySet());
	for(int i = 0; i < keys_list.size(); i++){

		String pheno = keys_list.get(i);
		Double[] input = new Double[6];

		LinkedList<Double>[] exist = pheno_to_all_pwy.get(pheno);

		for(int j = 0; j < exist.length; j++){

			double[] values = ArrayUtils.toPrimitive(exist[j].toArray(new Double[0]));

			if(values.length == 1){

				input[j] = values[0];
			}else{


				DescriptiveStatistics desc = new DescriptiveStatistics(values);
				input[j]=new Double(desc.getMean());
			}
			
		}	
		pheno_to_avr_pwy.put(pheno,input);
				
	}//i
}//method
		
}// inner class omics association
	
	
static class PO_GO {


static String[] po_cluster = {"defense","metabolism","phloem","protein_sorting","RNA","secretion","signal","sorting","transport","transporter","vesicular","xylem"};
static String[][] po_cluster_po_num;

static LinkedList<String>[][] cluster_genes_po;
	
void get_cluster_genes_po() throws IOException{

	System.out.println("get_cluster_genes:" );
	//TAIR	locus:2008960	SEC22		PO:0000013	TAIR:Publication:501715286|PMID:15806101	IEP		S	AT1G11890	AT1G11890|SEC22|ATSEC22|SECRETION 22|F12F1.27|F12F1_27	protein	taxon:3702	20081209	TAIR		TAIR:locus:2008960


	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern digit_pattern = Pattern.compile("\\d");

	cluster_genes_po  = new LinkedList[po_cluster.length][];
	po_cluster_po_num = new String[po_cluster.length][];


	String path = "data/cluster_po_analysis/";

	for(int i = 0; i < po_cluster.length; i++){

		String file_path = path + "po_anatomy_gene_arabidopsis_tair.assoc.nucleus.only." +  po_cluster[i];
		//System.out.println(po_cluster[i]);
		List<String> lines  = FileUtils.readLines(new File(file_path));
		LinkedList<String> po_num = new LinkedList();
		for(int q= 0; q < lines.size(); q++){
						
			//System.out.println(lines.get(q).trim());
			String line = lines.get(q).trim();
			String[] tab_split = tab_pattern.split(line);
			//String gene = tab_split[7].trim();
			String po = tab_split[3].trim();

			if(!po_num.contains(po)){

				po_num.add(po);
			}
			

		}//q

		po_cluster_po_num[i] = po_num.toArray(new String[0]);
		cluster_genes_po[i] = new LinkedList[po_cluster_po_num[i].length];

		for(int e = 0; e < po_cluster_po_num[i].length; e++){

			cluster_genes_po[i][e] = new LinkedList();

		}

		for(int q= 0; q < lines.size(); q++){
							
			String line = lines.get(q).trim();
			String[] tab_split = tab_pattern.split(line);
			String gene = tab_split[7].trim();
			String po = tab_split[3].trim();

			for(int e = 0; e < po_cluster_po_num[i].length; e++){

				if(po.equals(po_cluster_po_num[i][e])){

					if(!cluster_genes_po[i][e].contains(gene)){

						cluster_genes_po[i][e].add(gene);
						break;
					}
					
				}
			}
			
		}//q
				
	}//i
}

static String[] slim_cluster = {"abiotic","RNA","biotic","secretion","cytoskeleton","signal","Defense","skeleton","hormone","sorting","membrane_transport","transport","metabolism","vesicular","phloem","xylem","polarity"};

static String[][] slim_cluster_slim_num;

static LinkedList<String>[][] cluster_genes_slim;
	
void get_cluster_genes_slim() throws IOException{

	System.out.println("get_cluster_genes:" );
	//AT1G03950 locus:2024107 AT1G03950 involved in endosome transport via multivesicular body sorting pathway GO:0032509 27924 P other cellular processes IBA none PANTHER:PTN000049603|SGD:S000006435|TAIR:locus:2032607|PomBase:SPAC9E9.14|PomBase:SPAC4F8.01|UniProtKB:Q5AQN4|UniProtKB:Q9Y3E7|FB:FBgn0039402|TAIR:locus:2007883 Communication:501741973 2023-01-14



	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern digit_pattern = Pattern.compile("\\d");

	cluster_genes_slim  = new LinkedList[slim_cluster.length][];
	slim_cluster_slim_num = new String[slim_cluster.length][];


	String path = "data/cluster_slim_analysis/";

	for(int i = 0; i < slim_cluster.length; i++){

		String file_path = path + "ATH_GO_GOSLIM.nucleus." +  slim_cluster[i];
		
		List<String> lines  = FileUtils.readLines(new File(file_path));
		LinkedList<String> slim_num = new LinkedList();
		for(int q= 0; q < lines.size(); q++){
							
			String line = lines.get(q).trim();
			String[] space_split = space_pattern.split(line);
			//String gene = space_split[0].trim();
			String slim = "";
			for(int r = 0; r < space_split.length; r++){

				if(space_split[r].trim().startsWith("GO:")){

					slim = space_split[r].trim();
					break;
				}
			}

			if(!slim_num.contains(slim)){

				slim_num.add(slim);
			}
			

		}//q

		slim_cluster_slim_num[i] = slim_num.toArray(new String[0]);
		cluster_genes_slim[i] = new LinkedList[slim_cluster_slim_num[i].length];

		for(int e = 0; e < slim_cluster_slim_num[i].length; e++){

			cluster_genes_slim[i][e] = new LinkedList();

		}

		for(int q= 0; q < lines.size(); q++){
							
			String line = lines.get(q).trim();
			String[] space_split = space_pattern.split(line);
			String gene = space_split[0].trim();
			
			String slim = "";
			for(int r = 0; r < space_split.length; r++){

				if(space_split[r].trim().startsWith("GO:")){

					slim = space_split[r].trim();
					for(int e = 0; e < slim_cluster_slim_num[i].length; e++){

						if(slim.equals(slim_cluster_slim_num[i][e])){

							if(!cluster_genes_slim[i][e].contains(gene)){

								cluster_genes_slim[i][e].add(gene);
								break;
							}
							
						}
					}//e
					break;
				}
			}//r

		}//q
				
	}//i
}

 
static String[] go_tissue = {"flower","root","leaf","seed"};
static String[] go_score = {"queryitem",	"querytotal",	"bgitem",	"bgtotal",	"pvalue",	"FDR"};

static LinkedHashMap<String, String>[] go_term_type_to_term;
static LinkedHashMap<String, String>[] go_to_term_type;
static LinkedHashMap<String, Double[]>[] go_to_score;
static LinkedHashMap<String, LinkedList<String>>[] go_to_gene;
static LinkedHashMap<String, LinkedList<String>>[] gene_to_go;


void map_go_term_type_to_term() throws IOException{

//

	Pattern start_pattern = Pattern.compile("\\/\\/");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern semicolon_pattern = Pattern.compile(";");

	go_term_type_to_term = new LinkedHashMap[go_tissue.length];
	go_to_term_type = new LinkedHashMap[go_tissue.length];
	go_to_score = new LinkedHashMap[go_tissue.length];
	go_to_gene = new LinkedHashMap[go_tissue.length];
	gene_to_go = new LinkedHashMap[go_tissue.length];

	for(int i = 0; i < go_tissue.length; i++){
		go_term_type_to_term[i] = new LinkedHashMap();
		go_to_term_type[i] = new LinkedHashMap();
		go_to_score[i] = new LinkedHashMap();
		go_to_gene[i] = new LinkedHashMap();
		gene_to_go[i] = new LinkedHashMap();
	}

	String dir_path = "data/";

	for(int i = 0; i < go_tissue.length; i++){


		List<String> lines  = FileUtils.readLines(new File(dir_path + go_tissue[i] + ".GO.txt"));

		for(int q= 0; q < lines.size(); q++){

			String one_line = lines.get(q).trim();

			String[] tab_split = tab_pattern.split(one_line);
			String go_id = tab_split[0].trim();
			String term_type = tab_split[1].trim();
			String term = tab_split[2].trim();

			Double[] value = new Double[6];
			value[0] = new Double(tab_split[3].trim());
			value[1] = new Double(tab_split[4].trim());
			value[2] = new Double(tab_split[5].trim());
			value[3] = new Double(tab_split[6].trim());
			value[4] = new Double(tab_split[7].trim());
			value[5] = new Double(tab_split[8].trim());

			LinkedList<String> entries = new LinkedList();
			String[] entry_split = start_pattern.split(tab_split[9].trim());

			for(int w = 0; w < entry_split.length; w++){

				if(!entry_split[w].trim().equals("")){

					entries.add(entry_split[w].trim());
				}
			}
	
			go_term_type_to_term[i].put(term_type,term);
			go_to_term_type[i].put(go_id,term_type);
			go_to_score[i].put(go_id,value);
			go_to_gene[i].put(go_id,entries);
			
		}//q
	
		LinkedList<String> keys = new LinkedList<String>(go_to_gene[i].keySet());

		for(int z=0; z< keys.size(); z++){

			String go_id = keys.get(z);
			LinkedList<String> genes = go_to_gene[i].get(go_id);

			for(int q = 0; q < genes.size(); q++){

				String gene_id = genes.get(q).trim();

				if(gene_to_go[i].containsKey(gene_id)){

					LinkedList<String> exist = gene_to_go[i].get(gene_id);
					if(!exist.contains(go_id)){
						exist.add(go_id);
					}
					gene_to_go[i].put(gene_id,exist);
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(go_id);
					gene_to_go[i].put(gene_id,exist);
				}


	
			}//q

		}//z
	}//i
	
}//method


static LinkedHashMap<String,LinkedList<String>> atxg_to_po_ana;

void map_atxg_to_po_ana() throws IOException{
	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern at_pattern = Pattern.compile("@");
	Pattern comma_pattern = Pattern.compile(",");
	atxg_to_po_ana = new LinkedHashMap();

	
	String path = "data/out_po_get_po_atxg";

    	
		
    List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){
		//AT1G01010@[PO:0009005, PO:0009009, PO:0009029, PO:0009031, PO:0009047, PO:0009052, PO:0025022, PO:0000293]
	
		String one_line = lines.get(q).trim();
		String[] split_str = at_pattern.split(one_line);

		String atxg = split_str[0].trim();
		String po = split_str[1].trim();
		po = po.substring(1,po.length()-1);
		String[] comma_split = comma_pattern.split(po);

		LinkedList<String> po_list = new LinkedList();

		for(int e = 0; e < comma_split.length; e++){

			po_list.add(comma_split[e].trim());
		}

		atxg_to_po_ana.put(atxg,po_list);
		
			
 	}//q
 	
 }	

static LinkedHashMap<String,LinkedList<String>> atxg_to_po_tempo;

void map_atxg_to_po_tempo() throws IOException{
	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern at_pattern = Pattern.compile("@");
	Pattern comma_pattern = Pattern.compile(",");
	atxg_to_po_tempo = new LinkedHashMap();

	
	String path = "data/out_po_get_tempo_atxg";

    	
		
    List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){
		//AT1G01010@[PO:0009005, PO:0009009, PO:0009029, PO:0009031, PO:0009047, PO:0009052, PO:0025022, PO:0000293]
	
		String one_line = lines.get(q).trim();
		String[] split_str = at_pattern.split(one_line);

		String atxg = split_str[0].trim();
		String po = split_str[1].trim();
		po = po.substring(1,po.length()-1);
		String[] comma_split = comma_pattern.split(po);

		LinkedList<String> po_list = new LinkedList();

		for(int e = 0; e < comma_split.length; e++){

			po_list.add(comma_split[e].trim());
		}

		atxg_to_po_tempo.put(atxg,po_list);
		
			
 	}//q
 	
 }	

static LinkedHashMap<String,LinkedList<String>> atxg_to_po_slim;

void map_atxg_to_po_slim() throws IOException{
	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern at_pattern = Pattern.compile("@");
	Pattern comma_pattern = Pattern.compile(",");
	atxg_to_po_slim = new LinkedHashMap();

	
	String path = "data/out_get_slim_atxg";

    	
		
    List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){
		//AT1G01010@[GO:0005634, GO:0098542, GO:0043436, GO:0006355, GO:0009737, GO:0009617, GO:0010035, GO:0033993, GO:0006979, GO:0009414, GO:0003700, GO:0000976]
	
		String one_line = lines.get(q).trim();
		String[] split_str = at_pattern.split(one_line);

		String atxg = split_str[0].trim();
		String po = split_str[1].trim();
		po = po.substring(1,po.length()-1);
		String[] comma_split = comma_pattern.split(po);

		LinkedList<String> po_list = new LinkedList();

		for(int e = 0; e < comma_split.length; e++){

			po_list.add(comma_split[e].trim());
		}

		atxg_to_po_slim.put(atxg,po_list);
		
			
 	}//q
 	
 }//method
 
static LinkedList<String> homolog_list;
void get_homolog_gene() throws IOException{

	homolog_list = new LinkedList();

	
	String path = "data/araport.gene.geneName.homolog.sort.uniq3";
	List<String> all_lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < all_lines.size(); q++){

  		String one_line = all_lines.get(q).trim();
		homolog_list.add(one_line);
		
	}
}

static LinkedHashMap<String,String> gene_to_symbol;


void map_gene_to_symbol() throws IOException{

	gene_to_symbol = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/araport.gene.geneName.sort";
	List<String> all_lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < all_lines.size(); q++){

  		String one_line = all_lines.get(q).trim();
		String[] split_str = space_pattern.split(one_line);
		gene_to_symbol.put(split_str[0].trim(),split_str[1].trim());
		
	}
}//method

static LinkedHashMap<String,String> gene_to_ipd3;


void map_gene_to_ipd3() throws IOException{

	gene_to_ipd3 = new LinkedHashMap();
	Pattern space_pattern = Pattern.compile("\\s+");

	String path = "data/media-5.csv.gene.module.sort2";
	List<String> all_lines  = FileUtils.readLines(new File(path));
	for(int q= 0; q < all_lines.size(); q++){

  		String one_line = all_lines.get(q).trim();
		String[] split_str = space_pattern.split(one_line);
		gene_to_ipd3.put(split_str[0].trim(),split_str[1].trim());
		
	}
}//method

}//inner class PO_GO


static class medicago_cog_module {

static LinkedHashMap<String,String> medi_string_uniprot_to_gff_gene;
static LinkedHashMap<String,String> medi_gff_gene_to_string_uniprot;
static LinkedHashMap<String,String> medi_string_uniprot_to_gff_transcript;

void map_medi_string_uniprot_to_gff_transcript() throws IOException{

	medi_string_uniprot_to_gff_transcript = new LinkedHashMap();
	medi_string_uniprot_to_gff_gene = new LinkedHashMap();
	medi_gff_gene_to_string_uniprot = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("\\|");

	String path = "big/idmapping_2023_09_15.txt";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		if(one_line.contains("ID   ")){

			//System.out.println("498 one_line:" + one_line);
			String uni_id = "";
			String gene_id = "";
			String trans_id = "";
			String gn = "";

			for(int r = q+1; r < lines.size(); r++){
			
				String next_line = lines.get(r).trim();
				
				
				if(next_line.contains("ID   ")){
				
					break;
				}
			
				if(next_line.contains("AC   ")){
				
					uni_id = next_line.replaceFirst("AC   ","").trim();
					uni_id = uni_id.substring(0, uni_id.indexOf(";"));
				}else if(next_line.contains("GN   ")){
				
					
					gn = gn + next_line.substring(next_line.indexOf("   ")+1,next_line.length());
					
					
				}else if(next_line.contains("OS   ")){
				
					String[] str_split = space_pattern.split(gn.trim());
					
					for(int f =0; f < str_split.length; f++){
					
						if(str_split[f].trim().contains("=MTR_")){
						
							gene_id = str_split[f].trim().substring(str_split[f].trim().indexOf("=") +1,str_split[f].trim().length());
							break;
						}
					}
					
					//System.out.println("str_split.length:" + str_split.length);
					for(int f =0; f < str_split.length; f++){
					
						if(str_split[f].trim().contains("EMBL:") || str_split[f].trim().contains("EnsemblPlants:")){
						
							
							String[] bar_split = bar_pattern.split(str_split[f].trim());
							//System.out.println("EMBL:" + bar_split.length);
							
							for(int e = 0; e < bar_split.length; e++){
							
								if(bar_split[e].trim().contains("EMBL:") || bar_split[e].trim().contains("EnsemblPlants:")){
								
									//System.out.println("bar_split[e].trim():" + bar_split[e].trim());
									
									if(bar_split[e].trim().contains(".")){
										trans_id = bar_split[e].trim().substring(bar_split[e].trim().indexOf(":")+1,bar_split[e].trim().indexOf("."));
									}else{
									
										trans_id = bar_split[e].trim().substring(bar_split[e].trim().indexOf(":")+1,bar_split[e].trim().indexOf("}"));
									}
									medi_string_uniprot_to_gff_transcript.put(uni_id, trans_id);
									//System.out.println("error2:" + uni_id + ":" +  trans_id);
									break;
								}
							}
							
							break;
						}
					}
					
					
					
					
					medi_string_uniprot_to_gff_gene.put(uni_id, gene_id);
					//System.out.println("error:" + uni_id + ":" + gene_id);
					medi_gff_gene_to_string_uniprot.put(gene_id,uni_id);
					
					
					break;
				}

			}
		}
	}//lines

				
}//method

static LinkedHashMap<String,String> ara_string_uniprot_to_gff_transcript;

void map_ara_string_uniprot_to_gff_transcript() throws IOException{

	ara_string_uniprot_to_gff_transcript = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("\\|");

	String path = "data/Tair_atxg_nm_uniprot_id_map.idmaps.all";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		String[] space_split = space_pattern.split(one_line);
  		String uni_id = space_split[1].trim();
  		String mn_id = space_split[0].trim();
  		ara_string_uniprot_to_gff_transcript.put(uni_id,mn_id);
			
	}//lines

				
}//method


static LinkedHashMap<String,String> ara_string_uniprot_to_gff_gene;
static LinkedHashMap<String,String> ara_gff_gene_to_string_uniprot;

void map_ara_string_uniprot_to_gff_id() throws IOException{

	ara_string_uniprot_to_gff_gene = new LinkedHashMap();
	ara_gff_gene_to_string_uniprot = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("\\|");

	String path = "data/idmapping_2023_09_26.tsv.id.map";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
  		String[] space_split = space_pattern.split(one_line);
  		String uni_id = space_split[1].trim();
  		String atxg_id = space_split[0].trim();
		ara_string_uniprot_to_gff_gene.put(uni_id,atxg_id);
		ara_gff_gene_to_string_uniprot.put(atxg_id,uni_id);
		
			
	}//lines

				
}//method

static LinkedHashMap<String,String> ara_to_cog;
static LinkedHashMap<String,String> ara_to_cog_desc;
static LinkedHashMap<String,String> ara_cog_to_cog_desc;
static LinkedHashMap<String,LinkedList<String>> desc_to_ara_id;
static LinkedHashMap<String,String> desc_to_ara_cog_id;
static LinkedHashMap<String,LinkedList<String>> ara_cog_to_id;
static LinkedHashMap<String,LinkedList<String>> ara_cog_to_id_gff;

void map_ara_to_cog() throws IOException{

	ara_to_cog = new LinkedHashMap();
	ara_to_cog_desc = new LinkedHashMap();
	ara_cog_to_cog_desc = new LinkedHashMap();
	ara_cog_to_id = new LinkedHashMap();
	desc_to_ara_id = new LinkedHashMap();
	desc_to_ara_cog_id = new LinkedHashMap();
	ara_cog_to_id_gff = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("|");

	String path = "big/COG.mappings.v12.0.txt.arabidopsis";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
		//3702.Q8W498	1	1048	KOG2171	ARM repeat superfamily protein.
  		String[] str_split = tab_pattern.split(one_line);
  		String id = str_split[0].trim();
  		id = id.replaceFirst("3702.","");
  		
  		String ara_gff = ara_string_uniprot_to_gff_gene.get(id);
			
  		String cog = str_split[3].trim();
  		String desc = str_split[4].trim();
		ara_to_cog.put(id,cog);
		ara_to_cog_desc.put(id,desc);
		ara_cog_to_cog_desc.put(cog,desc);
		
		desc_to_ara_cog_id.put(desc,cog);
		
		if(desc_to_ara_id.containsKey(desc)){

			LinkedList<String> exist = desc_to_ara_id.get(desc);
			
			if(!exist.contains(id)){
			
				exist.add(id);
			}
 			desc_to_ara_id.put(desc,exist);
		}else{
		
			LinkedList<String> exist = new LinkedList();
			exist.add(id);
			desc_to_ara_id.put(desc,exist);
								
		}
		
		if(ara_cog_to_id.containsKey(cog)){

			LinkedList<String> exist = ara_cog_to_id.get(cog);
			
			if(!exist.contains(id)){
			
				exist.add(id);
			}
 			ara_cog_to_id.put(cog,exist);
		}else{
		
			LinkedList<String> exist = new LinkedList();
			exist.add(id);
			ara_cog_to_id.put(cog,exist);
								
		}
		
		if(ara_gff != null){
		
			if(ara_cog_to_id_gff.containsKey(cog)){

				LinkedList<String> exist = ara_cog_to_id_gff.get(cog);
				
				if(!exist.contains(ara_gff)){
				
					exist.add(ara_gff);
				}
	 			ara_cog_to_id_gff.put(cog,exist);
			}else{
			
				LinkedList<String> exist = new LinkedList();
				exist.add(ara_gff);
				ara_cog_to_id_gff.put(cog,exist);
									
			}
		}
		
	}//lines
	
	LinkedList<String> keys = new LinkedList<String>(ara_to_cog_desc.keySet());
	for(int j = 0; j < keys.size(); j++){
	
		String atxg_cog = keys.get(j);
		String desc = ara_to_cog_desc.get(atxg_cog);
		
		if(desc != null){
			//System.out.println("ara_desc:" + atxg_cog + ":" + ara_to_cog.get(atxg_cog) + ":" + desc);
		}
	}
		

				
}//method


static LinkedHashMap<String,String> medi_to_cog;
static LinkedHashMap<String,String> medi_to_cog_desc;
static LinkedHashMap<String,String> medi_cog_to_cog_desc;
static LinkedHashMap<String,LinkedList<String>> medi_cog_to_id;
static LinkedHashMap<String,LinkedList<String>> medi_cog_to_id_gff;
static LinkedHashMap<String,String> desc_to_medi_cog_id;


void map_medi_to_cog() throws IOException{

	medi_to_cog = new LinkedHashMap();
	medi_to_cog_desc = new LinkedHashMap();
	medi_cog_to_cog_desc = new LinkedHashMap();
	medi_cog_to_id = new LinkedHashMap();
	medi_cog_to_id_gff = new LinkedHashMap();
	desc_to_medi_cog_id = new LinkedHashMap();
	
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("|");

	String path = "big/COG.mappings.v12.0.txt.medicago";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
		//3702.Q8W498	1	1048	KOG2171	ARM repeat superfamily protein.
  		String[] str_split = tab_pattern.split(one_line);
  		String id = str_split[0].trim();
  		id = id.replaceFirst("3880.","");  		
  		
  		String medi_gff = medi_string_uniprot_to_gff_gene.get(id);		
		
  		String cog = str_split[3].trim();
  		String desc = str_split[4].trim();
  		medi_to_cog.put(id,cog);
		medi_to_cog_desc.put(id,desc);
		medi_cog_to_cog_desc.put(cog,desc);
		desc_to_medi_cog_id.put(desc,cog);
		
		if(medi_cog_to_id.containsKey(cog)){

			LinkedList<String> exist = medi_cog_to_id.get(cog);
			
			if(!exist.contains(id)){
			
				exist.add(id);
			}
 			medi_cog_to_id.put(cog,exist);
		}else{
		
			LinkedList<String> exist = new LinkedList();
			exist.add(id);
			medi_cog_to_id.put(cog,exist);
								
		}
  		
		if(medi_gff != null){
		
			if(medi_cog_to_id_gff.containsKey(cog)){

				LinkedList<String> exist = medi_cog_to_id_gff.get(cog);
				
				if(!exist.contains(medi_gff)){
				
					exist.add(medi_gff);
				}
	 			medi_cog_to_id_gff.put(cog,exist);
			}else{
			
				LinkedList<String> exist = new LinkedList();
				exist.add(medi_gff);
				medi_cog_to_id_gff.put(cog,exist);
									
			}
		}
		
	}//lines
	
	LinkedList<String> keys = new LinkedList<String>(medi_to_cog_desc.keySet());
	for(int j = 0; j < keys.size(); j++){
	
		String medi_cog = keys.get(j);
		String desc = medi_to_cog_desc.get(medi_cog);
		
		if(desc != null){
			//System.out.println("medi_desc:" + medi_cog + ":" + medi_to_cog.get(medi_cog) + ":" + desc);
		}
	}
				
}//method



static List<String> cog_ara_medi_intersection;

void get_cog_ara_medi_intersection() throws IOException{

	String path = "data/COG.mappings.v12.0.txt.arabidopsis.medicago.cod.id.intersection";

	cog_ara_medi_intersection  = FileUtils.readLines(new File(path));
}

static LinkedHashMap<LinkedList<String>,LinkedList<String>> match_ara_gff_to_medi_gff;

void match_cogs_between_ara_medi_intersection() throws IOException{

	match_ara_gff_to_medi_gff = new LinkedHashMap();
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("|");
	
	for(int q= 0; q < cog_ara_medi_intersection.size(); q++){
	
		String cog = cog_ara_medi_intersection.get(q);
		LinkedList<String> medi_gff = medi_cog_to_id_gff.get(cog);
		LinkedList<String> ara_gff = ara_cog_to_id_gff.get(cog);
		
		if(medi_gff != null && ara_gff != null){
		
			match_ara_gff_to_medi_gff.put(ara_gff,medi_gff);
			System.out.println("intersection cog:" + cog + ":" +   medi_cog_to_cog_desc.get(cog) + ":" + ara_cog_to_cog_desc.get(cog) + ":" + Arrays.toString(medi_gff.toArray()) + ":" + Arrays.toString(ara_gff.toArray()));
		}
	}
}



static LinkedList<String> ara_uniprot_module_e_cog;
static LinkedList<String> ara_uniprot_module_e_cog_gff;

void get_ara_uniprot_module_e_cog() throws IOException{

	ara_uniprot_module_e_cog = new LinkedList();
	ara_uniprot_module_e_cog_gff = new LinkedList();
	
	String path = "data/out_demonstration_medi_ara_cog_ipd3_metal.ara_medi_cog.module.e.sort.uniq.sel.ara.only.sort.uniq";

	List<String> lines  = FileUtils.readLines(new File(path));

	
	for(int q= 0; q < lines.size(); q++){

		String uni_id = lines.get(q).trim();
		ara_uniprot_module_e_cog.add(uni_id);
		
		String gff_id_temp = ara_string_uniprot_to_gff_gene.get(uni_id);
		
		if(gff_id_temp != null){
		
			ara_uniprot_module_e_cog_gff.add(gff_id_temp);
		}
	}

}


static LinkedList<String> medi_uniprot_module_e_cog;
static LinkedList<String> medi_gff_module_e_cog;

void get_module_e_medi_cog() throws IOException{

	medi_uniprot_module_e_cog = new LinkedList();
	medi_gff_module_e_cog = new LinkedList();
	
	String path = "data/out_demonstration_medi_ara_cog_ipd3_metal.ara_medi_cog.module.e.sort.uniq.sel.medi.only.sort.uniq";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String uni_id = lines.get(q).trim();
		medi_uniprot_module_e_cog.add(uni_id);
		String gff_id = medi_string_uniprot_to_gff_transcript.get(uni_id);
		
		if(gff_id != null){
		
			medi_gff_module_e_cog.add(gff_id);
		}
	}
		

}

static LinkedList<String> medi_deg_gene_id;
static LinkedList<String> medi_deg_uni_id;

void get_medi_deg() throws IOException{ 

	medi_deg_gene_id = new LinkedList();
	medi_deg_uni_id = new LinkedList();

	String path = "data/deg.medi.id";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String gene_id = lines.get(q).trim();
		medi_deg_gene_id.add(gene_id);
		String uni_id = medi_gff_gene_to_string_uniprot.get(gene_id);
		FileUtils.writeStringToFile(new File("data/deg.medi.id.uniprot"),uni_id+"\n",true);
		medi_deg_uni_id.add(uni_id);
	}
}

static LinkedList<String>[] medi_to_inter_part;

void map_medi_to_inter_part() throws IOException{

	medi_to_inter_part = new LinkedList[medi_deg_uni_id.size()];
	
	for(int i = 0; i < medi_to_inter_part.length; i++){
	
		medi_to_inter_part[i] = new LinkedList();
	}
	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("	");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	Pattern bar_pattern = Pattern.compile("|");

	String path = "big/3880.protein.links.v12.0.DEG.txt";

	List<String> lines  = FileUtils.readLines(new File(path));
	
	
	for(int q= 0; q < lines.size(); q++){

		String one_line = lines.get(q).trim();
		String[] str_split = space_pattern.split(one_line);
		
  		String id = str_split[0].trim();
  		int id_ind = -1;
  		
  		for(int z = 0; z < medi_deg_uni_id.size(); z++){
  		
  			if(id.equals(medi_deg_uni_id.get(z))){
  			
  				id_ind =z;
  				break;
  			}
  		}
  		
  		if(id_ind != -1){
  		
  			medi_to_inter_part[id_ind].add(str_split[1].trim());
  		}else{
  		
  			System.out.println("error interaction partner not found:" + one_line);
  		}
		
	}//lines
				
}//method

	
}//inner class medicago_cog_module


}//class
	
