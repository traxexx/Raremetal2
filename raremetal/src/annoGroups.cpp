#include "annoGroups.h"
#include "MetaUtility.h"
#include "SummaryFileReader.h"

// load group file to vectors
// name: annoGroups
// below vector<>
// chr: annoChrs;
// positions: annoPositions
// refs: annoRefs;
// alts: annoAlts;
void annoGroups::LoadGroupFile( String& group_file_name )
{
	cond_number = 0;
	Covs = NULL;
	// fill in in anno*
	FILE * file = fopen(group_file_name,"r");
	if(file==NULL)
		error("ERROR! Cannot open group file %s.\n",group_file_name.c_str());
	int idx = 0;
	while (!feof(file)) {
		StringArray tmp;
		String buffer;
		buffer.ReadLine(file);
		tmp.AddTokens(buffer, "\t ");
		if(tmp.Length()<=1)
			continue;
		annoGenes.push_back(tmp[0]);
		bool status = addMarkerGroupToAnno( 1,tmp );
		if (!status)
			error("Abnormal marker in group file:\n  %s",buffer.c_str());
		idx++;
	}
	fclose(file);
	// additional steps
	initializeMarkerIndex();
	initializeDataStorage();
}

void annoGroups::initializeDataStorage()
{
	groupUs = new Vector [annoGenes.size()];
	groupVs = new Vector [annoGenes.size()];
	Mafs = new Vector[annoGenes.size()];
	for(int g=0;g<annoGenes.size();g++) {
		groupUs[g].Dimension(annoChrs[g].size(),0);
		groupVs[g].Dimension(annoChrs[g].size(),0);
		Mafs[g].Dimension(annoChrs[g].size(),0);
	}
}

void annoGroups::initializeMarkerIndex()
{
	if (annoChrs.empty())
		error("No valid group information!\n");
	for(int i=0;i<annoChrs.size();i++) {
		for(int j=0;j<annoChrs[i].size();j++) {
			if (markerIndex.find(annoChrs[i][j])==markerIndex.end())
				markerIndex[annoChrs[i][j]];
			String key = "";
			key += annoPositions[i][j];
			if (newFormat)
				key += ":" +  annoRefs[i][j] + ":" + annoAlts[i][j];
			std::pair<int,int> tmp (-1,-1);
			markerIndex[annoChrs[i][j]][key] = tmp;
		}
	}
}

bool annoGroups::addMarkerGroupToAnno( int offset, StringArray& tmp )
{
	std::vector<String> chrs;
	std::vector<int> positions;
	std::vector<String> refs;
	std::vector<String> alts;
	for(int i=offset;i<tmp.Length();i++) {
		StringArray tmp2;
		tmp2.AddTokens(tmp[i],":");
		if (tmp2.Length() != 4)
			return false;
		chrs.push_back(tmp2[0]);
		positions.push_back(tmp2[1]);
		refs.push_back(tmp2[2]);
		alts.push_back(tmp2[3]);
	}
	annoChrs.push_back(chrs);
	annoPositions.push_back(positions);
	annoRefs.push_back(refs);
	annoAlts.push_back(alts);
	return true;
}


// load Cov Files
// store in markers***
// only load marker line that exist in group file
void annoGroups::LoadCovFile(String& cov_file_name,int _adjust, int _marker_col, int _cov_col, int _sampleSize)
{
	Covs = NULL;
	adjust = _adjust;
	marker_col = _marker_col;
	cov_col = _cov_col;
	sampleSize = _sampleSize;
	loadCovStrings(cov_file_name);
	// now set marker positions
	setMarkersInCovs();
	// clear & additional steps
	clearLoadData();
}

void annoGroups::loadCovStrings(String& cov_file_name)
{
	// check if clear before start
	if (!markersExp.empty())
		error("Data structure not cleared!\n\n");

	// start	
	SummaryFileReader covReader;
	printf("Reading cov matrix from %s ...\n",cov_file_name.c_str());
	IFILE covfile_;
	covfile_  = ifopen(cov_file_name,"r");
	if(covfile_ == NULL)
		error("ERROR! Cannot open file: %s! Input cov file has to be bgzipped and tabix indexed using the following command:\n bgzip yourfile.singlevar.cov.txt; tabix -c \"#\" -s 1 -b 2 -e 2 yourprefix.singlevar.cov.txt.gz\n",cov_file_name.c_str());
//	Tabix covtabix;
//	String tabix_name = filename + ".tbi";
//	StatGenStatus::Status libstatus = covtabix.readIndex( tabix_name.c_str() );

	int m=0;
	bool pass_header = 0;
	int vec_idx = 0;
	int genome_idx = 0;
	while (!ifeof(covfile_)) {
		String buffer;
		buffer.ReadLine(covfile_);
		if (!pass_header) {
			if (buffer.Find("CHROM")==-1)
				continue;
			// now check new or old format
			StringArray tokens;
			tokens.AddTokens(buffer,"\t ");
			if (tokens[2]=="MARKERS_IN_WINDOW" && tokens[3]=="COV_MATRICES")
				newFormat = false;
			else if (tokens[2]=="REF" && tokens[3]=="ALT" && tokens[4]=="EXP" && tokens[5]=="COV_MATRICES")
				newFormat = true;
			else
				error("Header line: %s\nCovariance matrix is neither new or old format...are you using the right file?\n\n",buffer.c_str());
			pass_header = 1;
			continue;
		}
		StringArray tokens;
		tokens.AddTokens(buffer,"\t ");
//		if ( RegionStatus ) {
//			if ( tokens[1].AsInteger() > End || tokens[0] != Chr ) // out of this region or into another chromosome
//				break;
//		}
		m++;
		bool flip_status = false;
		int position_int = tokens[1].AsInteger();
		bool in_status = isInMarkerIndex(tokens, flip_status);
		if (!in_status) { // variant does not show up in group file or cond analysis
			genome_idx++;
			continue;
		}
		// now add
		removeChrFromString(tokens[0]);
		String key = tokens[1];
		if (newFormat)
			key += ":" + tokens[2] + ":" + tokens[3];
		markerIndex[tokens[0]][key].first = genome_idx;
		markerIndex[tokens[0]][key].second = vec_idx;
		if (newFormat) {
			int exp_col = 4;
			int new_cov_col = cov_col+2;
			markersExp.push_back(tokens[exp_col]);
			markersCov.push_back(tokens[new_cov_col]);
		}
		else { // additionally check flip_allele_map
			String key = tokens[0] + ":" + tokens[1];
			if (flip_allele_map.find(key)!=flip_allele_map.end())
				flip_status = true;
			markersCov.push_back(tokens[cov_col]);
		}
		markersFlip.push_back(flip_status);
		genome_idx++;
		vec_idx++;
	}
	ifclose(covfile_);
	printf("done\n");		
}


bool annoGroups::isInMarkerIndex( StringArray& tokens, bool& flip_status )
{
	flip_status = false;
	if (markerIndex.find(tokens[0])==markerIndex.end())
		return false;
	String key;
	if (newFormat)
		key = tokens[1] + ":" + tokens[2] + ":" + tokens[3];
	else
		key = tokens[1];
	if (newFormat) {
		if (markerIndex[tokens[0]].find(key)==markerIndex[tokens[0]].end())
			return true;
		else {
			key = tokens[1] + ":" + tokens[3] + ":" + tokens[2]; 
			if (markerIndex[tokens[0]].find(key)!=markerIndex[tokens[0]].end()) {
				flip_status = true;
				return true;
			}
			else
				return false;
		}
	}
	else {
		if (markerIndex[tokens[0]].find(key)!=markerIndex[tokens[0]].end())
			return true;
		else
			return false;
	}
}

void annoGroups::setMarkersInCovs()
{
	// initialize
	if (Covs==NULL)
		Covs = new Matrix [annoGenes.size()];

	// extract
	for(int g=0;g<annoGenes.size();g++) {
		if (groupUs[g].Length()==0)
			printf("Group %s has no variant. Do not set marker covs!\n",annoGenes[g].c_str());
		int n = cond_number>0 ? cond_number : annoPositions[g].size();
		Covs[g].Dimension(n,n,0);
		for(int i=0;i<n;i++) {
			StringArray covs;
			int pi = annoPositions[g][i];
			String name_i = "";
			name_i += pi;
			if (newFormat)
				name_i += ":" + annoRefs[g][i] + ":" + annoAlts[g][i];
			int i1 = markerIndex[annoChrs[g][i]][name_i].first;
			int i2 = markerIndex[annoChrs[g][i]][name_i].second;
			if (i1==-1&&i2==-1)
				continue;
			else if (i1==-1||i2==-1)
				error("At marker %s:%s,gen_idx=%d,vec_idx=%d. Something is wrong!\n\n",annoChrs[g][i].c_str(),name_i.c_str(),i1,i2);
			if (i2>=markersCov.size())
				error("At marker %s:%s,vec_idx=%d, which is larger than cov size of %d. Something is wrong!\n\n",annoChrs[g][i].c_str(),name_i.c_str(),i2,markersCov.size());
			if (newFormat)
				addNewFormatCov( markersExp[i2], markersCov[i2],covs);
			else
				covs.AddTokens( markersCov[i2],",");
			double flip_i = 1;
			if (markersFlip[i2])
				flip_i = -1;
			int nj = cond_number>0 ? marker_number : n;
			for(int j=i;j<nj;j++) {
				int pj = annoPositions[g][j];
				String name_j = "";
				name_j += pj;
				if (newFormat)
					name_j += ":" + annoRefs[g][j] + ":" + annoAlts[g][j];
				int j1 = markerIndex[annoChrs[g][i]][name_j].first;
				int j2 = markerIndex[annoChrs[g][i]][name_j].second;
				if (j1==-1&&j2==-1)
					continue;
				double flip_j = 1;
				if (markersFlip[j2])
					flip_j = -1;
				int idx = j1-i1;
				if (idx<0)
					idx *= (-1);
				double val;
				if (idx>covs.Length())
					val = 0;
				else {
					val = flip_i*flip_j * covs[idx].AsDouble()*sampleSize;
					if (newFormat>0)
						val *= pow(10,markersExp[i2]);
				}
				Covs[g][i][j] += val;
			}
		}
	}
}

void annoGroups::clearLoadData()
{
	// clear after filling up each study
	for(std::map< String, std::map<String, std::pair<int,int> > >::iterator p1=markerIndex.begin();p1!=markerIndex.end();p1++) {
		for(std::map<String, std::pair<int,int> >::iterator p2=p1->second.begin();p2!=p1->second.end();p2++) {
			p2->second.first = -1;
			p2->second.second = -1;
		}
	}
	markersExp.clear(); // for new format
	markersInWindow.Clear(); // for old format
	markersCov.clear();	
	markersFlip.clear();
}

// for new format
// read marker cov then add to vector
void annoGroups::addNewFormatCov( int mexp, String & cov_str, StringArray& covs)
{
	StringArray commas;
	commas.AddTokens(cov_str,',');
	// length of covs
	int n = commas.Length();
	if (n<1)
		error("At line: %s:...,no index of covariance matrices! Are you using the right cov file?\n",cov_str.c_str());
	bool index_exist = 0;
	// now check if index exists
	StringArray first_tokens;
	first_tokens.AddTokens(commas[0],":");
	if (first_tokens.Length()==2)
		index_exist = 1;
	else if (first_tokens.Length()!=1)
		error("At line: %s:...,abnormal index of covariance matrices! Are you using the right cov file?\n",cov_str.c_str());

	// if no index, add directly
	if (index_exist) {
		int cov_len = -1;
		for(int i=n-1;i>=0;i--) {
			StringArray tokens;
			tokens.AddTokens( commas[i],':' );
			if (tokens.Length()==1)
				continue;
			cov_len = tokens[1].AsInteger()+1;
			break;
		}
		if (cov_len==-1)
			error("At line: %s:...,no index of covariance matrices! Are you using the right cov file?\n",commas[0].c_str());
		covs.Dimension(cov_len);
		// now add
		int last_index = 0;
		for(int i=0;i<n;i++) {
			StringArray tokens;
			tokens.AddTokens( commas[i],':' );		
			if (tokens.Length()>2||tokens.Length()<1)
				error("At line: ...:%s:...,abnormal separation!\n",commas[i].c_str());
			if (tokens.Length()==1)
				last_index++;
			else // tokens ==2
				last_index = tokens[1].AsInteger();
	//printf("covlen=%d,last_index=%d\n",cov_len,last_index);		
			covs[last_index] = tokens[0];
		}
	}
	else {
		covs.Dimension(n);
		for(int i=0;i<n;i++)
			covs[i] = commas[i];
	}
}

// load conditional analysis file
// set conditional analysis markers as groups (same as group test)
void annoGroups::LoadCondFile( String& cond_file_name )
{
	IFILE condFile = ifopen(cond_file_name,"r");
	if(condFile==NULL)
		error("Can not open file %s.\n",cond_file_name.c_str());
	annoGenes.push_back("NA");
	annoChrs.resize(1);
	annoPositions.resize(1);
	annoRefs.resize(1);
	annoAlts.resize(1);
	int idx = 0;
	while (!ifeof(condFile)) {
		String buffer;
		buffer.ReadLine(condFile);
		if(buffer.FindChar('#')!=-1 || buffer.FindChar(':')==-1)
			continue;
		StringArray tmpMarker;
		tmpMarker.AddTokens(buffer, ":");
		if (tmpMarker.Length()!=4)
			error("Abnormal line at conditional file:\n%s",buffer.c_str());
		annoChrs[0].push_back(tmpMarker[0]);
		annoPositions[0].push_back(tmpMarker[1]);
		annoRefs[0].push_back(tmpMarker[2]);
		annoAlts[0].push_back(tmpMarker[3]);
		idx++;
	}
	ifclose(condFile);
	cond_number = idx;
	initializeMarkerIndex();
	initializeDataStorage();
}

void annoGroups::LoadSingleMarker(int g,String& chr,int position,String& ref,String& alt)
{
	if (annoGenes.size()<=g)
		error("[annoGroups::LoadSingleMarker] annoGenes size = %d but g = %d\n\n",annoGenes.size(),g);
	annoChrs[g].push_back(chr);
	annoPositions[g].push_back(position);
	annoRefs[g].push_back(ref);
	annoAlts[g].push_back(alt);
}

void annoGroups::FlipAllele(int g,int i)
{
	String tmp = annoRefs[g][i];
	annoRefs[g][i] = annoAlts[g][i];
	annoAlts[g][i] = tmp;
	String key = annoChrs[g][i] + ":";
	key += annoPositions[g][i];
	flip_allele_map[key] = true;
}

void annoGroups::SetU(int g, int i, double u)
{
	groupUs[g][i] = u;
}

void annoGroups::SetV(int g, int i, double v2)
{
	groupVs[g][i] = v2;
}

void annoGroups::SetMaf(int g, int i, double maf)
{
	Mafs[g][i] = maf;
}


String annoGroups::GetMarkerName(int g, int i)
{
	String marker = annoChrs[g][i] + ":" + annoPositions[g][i] + ":" + annoRefs[g][i] + ":" + annoAlts[g][i];
	return marker;
}

int annoGroups::VariantCountInGenes(int g)
{
	return annoChrs[g].size();
}

int annoGroups::GetGeneNumber()
{
	return annoGenes.size();
}

String annoGroups::GetChr(int g,int i)
{
	return annoChrs[g][i];
}

int annoGroups::GetPosition(int g, int i)
{
	return annoPositions[g][i];
}

String annoGroups::GetRef(int g,int i)
{
	return annoRefs[g][i];
}

String annoGroups::GetAlt(int g,int i)
{
	return annoAlts[g][i];
}

double annoGroups::GetOneCov(int g,int m1,int m2)
{
	return Covs[g][m1][m2];
}

void annoGroups::UpdateCovValue(double new_val,int g,int m1,int m2)
{
	Covs[g][m1][m2] = new_val;
}

// export matrix* cov of a single group
void annoGroups::ExportCov( Matrix& cov, int g )
{
	int n = Covs[g].rows;
	int m = Covs[g].cols;
	cov.Dimension(n,m);
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			cov[i][j] = Covs[g][i][j];
}


void annoGroups::MakeFullMatrix()
{
	// make half matrix as full matrix
	if (is_cov_half_matrix&&cond_number==0) {
		for(int g=0;g<annoGenes.size();g++) {
			for(int i=0;i<Covs[g].rows;i++)
				for(int j=i;j<Covs[g].rows;j++)
					 Covs[g][j][i] = Covs[g][i][j];
		}
	}
	is_cov_half_matrix = false;
}


// group test part
void annoGroups::GroupTest(String& method, String& prefix)
{
	runBurdenTest(method);
/*	if (method=="VT")
		runVt();
	else if (method=="SKAT")
		runSKAT();
	else { // regular burden test
		runBurdenTest(method);
	}*/
	String filename = prefix + ".meta." + method + ".results";
	IFILE f = ifopen(filename,"w");
	ifprintf(f,"#GROUP\tCOUNT\tVARIANTS\tAVR_MAF\tMIN_MAF\tMAX_MAF\tEFF_SIZE\tP_VALUE\n");
	for(int g=0;g<annoGenes.size();g++) {
		printGroupResult( g,f );	
	}
	ifclose(f);
	printf("  done.\n\n");
}

// remove the missing data in groupU, groupV and Cov
void annoGroups::RemoveMissingData()
{
	// check missingness of genes and single variants
	for(int g=0;g<annoGenes.size();g++) {
		std::vector<int> valid;
		for(int i=0;i<groupVs[g].Length();i++) {
			if (groupVs[g][i]!=0)
				valid.push_back(i);
		}
		if (valid.empty()) {
			groupVs[g].Dimension(0);
			continue;
		}
		if (valid.size()<groupVs[g].Length()) {
			// remove only singles
			for(int v=0;v<valid.size();v++) {
				int old_idx = valid[v];
				annoChrs[g][v] = annoChrs[g][old_idx];
				annoPositions[g][v] = annoPositions[g][old_idx];
				annoRefs[g][v] = annoRefs[g][old_idx];
				annoAlts[g][v] = annoAlts[g][old_idx];
				groupUs[g][v] = groupUs[g][old_idx];
				groupVs[g][v] = groupVs[g][old_idx];
				Mafs[g][v] = Mafs[g][old_idx];
				Matrix tmp = Covs[g];
				for(int a=0;a<valid.size();a++) // do it per row
					Covs[g][v][a] = tmp[old_idx][valid[a]];
			}
			annoChrs[g].resize(valid.size());
			annoPositions[g].resize(valid.size());
			annoRefs[g].resize(valid.size());
			annoAlts[g].resize(valid.size());
			groupUs[g].Dimension(valid.size());
			groupVs[g].Dimension(valid.size());
			Mafs[g].Dimension(valid.size());
			Covs[g].Dimension(valid.size(),valid.size());
		}
	}

	// search for empty elements
	std::vector<int> valid;
	for(int g=0;g<annoGenes.size();g++) {
		if (groupVs[g].Length()!=0)
			valid.push_back(g);
		else
			printf("Warning: Gene #%d %s with no available MAF is skipped for group test!\n",g+1,annoGenes[g].c_str());
	}
	if (valid.size()==annoGenes.size()) { // no missing element
		// additional initialization
		metaUs.Dimension(annoGenes.size(),0);
		metaVs.Dimension(annoGenes.size(),0);
		return;
	}
	// delete element if it is empty
	Matrix* p = Covs;
	Vector* pu = groupUs;
	Vector* pv = groupVs;
	Vector* pmaf = Mafs;
	groupUs = new Vector [valid.size()];
	groupVs = new Vector [valid.size()];
	Mafs = new Vector [valid.size()];
	Covs = new Matrix [valid.size()];
	for(int i=0;i<valid.size();i++) {
		int old_idx = valid[i];
		annoGenes[i] = annoGenes[old_idx];
		annoChrs[i] = annoChrs[old_idx];
		annoPositions[i] = annoPositions[old_idx];
		annoRefs[i] = annoRefs[old_idx];
		annoAlts[i] = annoAlts[old_idx];
		groupUs[i].Dimension(pu[old_idx].Length());
		groupUs[i] = pu[old_idx];
		groupVs[i].Dimension(pu[old_idx].Length());
		groupVs[i] = pv[old_idx];
		Mafs[i].Dimension(pu[old_idx].Length());
		Mafs[i] = pmaf[old_idx];
		Covs[i] = p[old_idx];
	}
	annoGenes.resize(valid.size());
	annoChrs.resize(valid.size());
	annoPositions.resize(valid.size());
	annoRefs.resize(valid.size());
	annoAlts.resize(valid.size());
	delete [] pu;
	delete [] pv;
	delete [] pmaf;
	delete [] p;

	// additional initialization
	metaUs.Dimension(annoGenes.size(),0);
	metaVs.Dimension(annoGenes.size(),0);
}

void annoGroups::setWeight(String& method, int g)
{
	Weight.Dimension(groupUs[g].Length());
	if(method=="burden")// equal weight
		for(int w=0;w<Weight.Length();w++)
				Weight[w] = 1.0;
	else if(method=="MB") // weight by 1/sqrt( maf*(1-maf) )
		for(int w=0;w<Weight.Length();w++)
			Weight[w] = sqrt(Mafs[g][w]*(1.0-Mafs[g][w]));
	else if (method=="MAB") {
		for(int w=0;w<Weight.Length();w++)
			Weight[w] = Mafs[g][w];
	}
	else if (method=="BBeta") { // truncated beta
		double alpha = 0.5;
		double beta = 0.5;
		double f0 = 2 / sampleSize; // truncate at 4 alleles
		for(int w=0;w<Weight.Length();w++) {
			double xmaf = Mafs[g][w];
			if (xmaf>0 && xmaf<f0)
				xmaf = f0;
			if (xmaf>0 && xmaf>(1-f0))
				xmaf = 1-f0;
			double beta_density = GetBetaDensity(alpha,beta,xmaf);
			Weight[w] = (beta_density*beta_density);
		}
	}
	else
		error("Invalid weight %s!\n",method.c_str());
	// for burden test, need to 1/w
	for(int w=0;w<Weight.Length();w++)
		Weight[w] = 1/Weight[w];
}

void annoGroups::runBurdenTest(String& method)
{
	for(int g=0;g<annoGenes.size();g++) {
		if (Mafs[g].Length()==0)
			error("At #%d gene, no available MAF. Something is wrong!\n",g+1);

		setWeight(method,g);
		
		double numerator  = Weight.InnerProduct(groupUs[g]);
		Vector tmp;
		tmp.Dimension(annoChrs[g].size());
		
		for(int i=0;i<tmp.Length();i++)
			tmp[i] = Weight.InnerProduct(Covs[g][i]);
		double denominator = tmp.InnerProduct(Weight);

		metaUs[g] = numerator;
		metaVs[g] = denominator;
	}		
}


// print one line of result
void annoGroups::printGroupResult(int g, IFILE& f)
{
	ifprintf(f, "%s\t%d\t", annoGenes[g].c_str(), annoChrs[g].size());
	for(int i=0;i<annoChrs[g].size();i++) {
		if (i>0)
			ifprintf(f,",");
		ifprintf(f,"%s:%d:%s:%s", annoChrs[g][i].c_str(), annoPositions[g][i],annoRefs[g][i].c_str(), annoAlts[g][i].c_str());
	}
	double average_af = Mafs[g].Average();
	double min_af = Mafs[g].Min();
	double max_af = Mafs[g].Max();
	ifprintf(f,"\t%g\t%g\t%g", average_af, min_af, max_af);

	if (metaVs[g]==0.0) {
		ifprintf(f,"\tNA\tNA\tNA\n");
		return;
	}

	double chisq = metaUs[g]*metaUs[g]/metaVs[g];
	double pvalue = pchisq(chisq,1,0,0);
	double effSize = metaUs[g]/metaVs[g];

	bool disect=false;
	while(pvalue==0.0) {
		disect=true;
		chisq *= 0.999;
		pvalue = pchisq(chisq,1,0,0);
	}

	ifprintf(f,"\t%g\t%g", effSize, disect?"<":"",pvalue);
	ifprintf(f,"\n");
}


/* variant threshold test
void annoGroups::runVt()
{
	printf("Performing Variable Threshold tests ...\n");
	printVTheader();
	
	StringArray chr_plot;
	Vector pos_plot;
	StringArray geneLabels;
	for(int g=0;g<annoGenes.size();g++) {
		if (Mafs[g].Length()==0)
			error("At #%d gene, no available MAF. Something is wrong!\n",g+1);
		if(g>1 && g%1000==1)
			printf("Finished analyzing %d genes.\n",g-1);
		
		Vector pvalue_VT;
		if(Mafs[g].Length()==1) {
			double chisq = groupUs[g][0]*groupUs[g][0]/groupVs[g][0];
			double p = pchisq(chisq,1,0,0);
			pvalue_VT.Push(p);
			ifprintf(output,"%s\t1\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),group.SNPlist[g][0].c_str(),maf[g][0],singleEffSize[g][0],singlePvalue[g][0],maf[g][0],maf[g][0],maf[g][0],singleEffSize[g][0],maf[g][0],singlePvalue[g][0]);
			ifprintf(output,"\n");
		}
		else // standard vt
			VTassocSingle(g);

		chr_plot.Push(annoChrs[g]);
		pos_plot.Push(annoPositions[g]);
		geneLabels.Push(annoGenes[g]);
	}
	
	String name = "VT (maf<";
	name +=  MAF_cutoff;
	name +=  ")";
	String extraname = "";
	String demo="";
	double GC = GetGenomicControlFromPvalue(pvalue_VT);
	demo="GC=";
	demo += GC;
	writepdf.Draw(pdf,geneLabels,pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);
	if(cond!="")
	{
		name += " Conditional Analysis";
		double GC = GetGenomicControlFromPvalue(cond_pvalue_VT);
		demo="GC=";
		demo += GC;
		writepdf.Draw(pdf,geneLabels,cond_pvalue_VT,chr_plot,pos_plot,name,extraname,demo,true);
	}
	
	ifclose(output);
	printf("Done.\n\n");
}


// vt for a single group
double Meta::VTassocSingle(int g)
{
	//STEP1: sort maf[g] to find # of cutoffs
	Vector cp_maf;
	cp_maf.Copy(Mafs[g]);
	cp_maf.Sort();
	Vector maf_cutoff;
	maf_cutoff.Push(cp_maf[0]);	
	for(int i=1;i<cp_maf.Length();i++)
	{
		if(cp_maf[i]>maf_cutoff[maf_cutoff.Length()-1])
			maf_cutoff.Push(cp_maf[i]);
	} //now unique maf cutoffs are saved in maf_cutoff.
	double pvalue=_NAN_,chosen_cutoff=_NAN_,chosen_effSize=_NAN_;
	double numerator=0.0,denominator=0.0,t_max=_NAN_;
	Vector weight,tmp,chosen_weight,score;
	Matrix cov_weight;
	weight.Dimension(maf[g].Length());
	tmp.Dimension(annoChrs[g].Length());
	cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
	score.Dimension(maf_cutoff.Length());

	// STEP2: loop through each cutoff to find max score
	for(int i=0;i<maf_cutoff.Length();i++) {
		for(int w=0;w<weight.Length();w++) {
			if(maf[g][w]<=maf_cutoff[i])
				weight[w]=1.0;
			else
				weight[w]=0.0;
			cov_weight[i][w]=weight[w];
		}
		numerator = weight.InnerProduct(groupUs[g]);
		for(int d=0;d<tmp.Length();d++)
				tmp[d] = weight.InnerProduct(Covs[g][d]);
		denominator = tmp.InnerProduct(weight);
		
		if(denominator != 0.0) {
			double t_stat = fabs(numerator/sqrt(denominator));
			score[i]=t_stat;
			if(t_max==_NAN_) {
				t_max = t_stat;
				chosen_cutoff = maf_cutoff[i];
				chosen_weight.Copy(weight);
				chosen_effSize = numerator/denominator;
			}
			else {
				if(t_stat>t_max) {
					t_max = t_stat;
					chosen_cutoff = maf_cutoff[i];
					chosen_weight.Copy(weight);
					chosen_effSize = numerator/denominator;
				}
			}
		}
		else
			score[i]=0.0;
	}
	if(score.Max()==0.0) {
		printf("Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
		fprintf(log,"Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
		return pvalue;
	}
	Vector tmp_maf,tmp_eff,tmp_pvalue;
	for(int i=0;i<maf[g].Length();i++) {
		if(chosen_weight[i]==1.0)
			tmp_maf.Push(maf[g][i]);
	}
	for(int i=0;i<maf[g].Length();i++) {
		if(chosen_weight[i]==1.0)
		{
			tmp_eff.Push(singleEffSize[g][i]);
			tmp_pvalue.Push(singlePvalue[g][i]);
		}
	}
	
	double average_af = tmp_maf.Average();
	double min_af = tmp_maf.Min();
	double max_af = tmp_maf.Max();
	String var;
	for(int i=0;i<maf[g].Length()-1;i++) {
		if(chosen_weight[i]==1.0)
			var += group.SNPlist[g][i] + ";";
	}
	if(chosen_weight[maf[g].Length()-1]==1.0)
		var += group.SNPlist[g][maf[g].Length()-1];
	
	//STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
	Matrix cov_U,cov_U_tmp;
	cov_U_tmp.Product(cov_weight,cov[g]);
	Matrix cov_weight_trans;
	cov_weight_trans.Transpose(cov_weight);
	cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
	//Calculate covariance matrix for (T_1 ... T_#cutoff)
	Matrix cov_T;
	cov_T.Dimension(cov_U.rows,cov_U.cols);
	cov2cor(cov_U,cov_T);
	//STEP4: calculate VT pvalue and report.
	int cutoff = maf_cutoff.Length();
	double * lower = new double [cutoff];
	double * upper = new double [cutoff];
	double * mean = new double [cutoff];
	
	for(int i=0;i<cutoff;i++) {
		mean[i] = 0.0;
		lower[i] = -t_max;
		upper[i] = t_max;
	}
	
	//Use pmvnorm to calculate the asymptotic p-value
	Vector result;
	pmvnorm(lower,upper,mean,cov_T,false,result);
	
	if(result[0]==-1.0)
	{
		if(!condition)
		{
			ifprintf(output,"\n");
		}
		else
		{
			ifprintf(output,"%g\t%g\tERROR:CORR_NOT_POS_SEMI_DEF\n",chosen_effSize,chosen_cutoff);
		}
	}
	else
	{
		if(1.0-result[0]==0.0)
		{
			//           printf("gene %s has result %g\n",group.annoGroups[g].c_str(),1.0-result[0]);
			printf("Using Shuang's algorithm to calculate MVN pvalue for gene %s ... ",group.annoGroups[g].c_str());
				if(maf_cutoff.Length()>20)
			{
				while(maf_cutoff.Length()>20)
					maf_cutoff.Delete(0);
				
				double numerator,denominator,t_max=_NAN_;
				Vector weight,tmp,chosen_weight,score;
				Matrix cov_weight;
				weight.Dimension(maf[g].Length());
				tmp.Dimension(group.SNPlist[g].Length());
				cov_weight.Dimension(maf_cutoff.Length(),maf[g].Length());
				for(int i=0;i<maf_cutoff.Length();i++)
				{
					for(int w=0;w<weight.Length();w++)
					{
						if(maf[g][w]<=maf_cutoff[i])
							weight[w]=1.0;
						else
							weight[w]=0.0;
						cov_weight[i][w]=weight[w];
					}
					if(condition)
						numerator = weight.InnerProduct(cond_stats[g]);
					else
						numerator = weight.InnerProduct(stats[g]);
					
					for(int d=0;d<tmp.Length();d++)
							tmp[d] = weight.InnerProduct(cov[g][d]);
					denominator = tmp.InnerProduct(weight);
					if(denominator != 0)
					{
						double t_stat = fabs(numerator/sqrt(denominator));
						score.Push(t_stat);
						if(t_max==_NAN_)
						{
							t_max = t_stat;
							chosen_cutoff = maf_cutoff[i];
							chosen_weight.Copy(weight);
							chosen_effSize = numerator/denominator;
						}
						else
						{
							if(t_stat>t_max)
							{
								t_max = t_stat;
								chosen_cutoff = maf_cutoff[i];
								chosen_weight.Copy(weight);
								chosen_effSize = numerator/denominator;
							}
						}
					}
					else
						score.Push(0.0);
				}
				if(score.Max()==0.0) {
					printf("Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
					fprintf(log,"Warning: group %s does not have qualified variants to group.\n",group.annoGroups[g].c_str());
					return pvalue;
					printf("completed!\n");
				}
				Vector tmp_maf,tmp_eff,tmp_pvalue;
				for(int i=0;i<maf[g].Length();i++)
				{
					if(chosen_weight[i]==1.0)
						tmp_maf.Push(maf[g][i]);
				}
				
				for(int i=0;i<maf[g].Length();i++)
				{
					if(chosen_weight[i]==1.0) {
						tmp_eff.Push(singleEffSize[g][i]);
						tmp_pvalue.Push(singlePvalue[g][i]);
					}
				}
				average_af = tmp_maf.Average();
				min_af = tmp_maf.Min();
				max_af = tmp_maf.Max();
				
				String var;
				for(int i=0;i<maf[g].Length()-1;i++) {
					if(chosen_weight[i]==1.0)
						var += group.SNPlist[g][i] + ";";
				}
				if(chosen_weight[maf[g].Length()-1]==1.0)
					var += group.SNPlist[g][maf[g].Length()-1];
				//STEP3: calculate covariance matrix for (U_1 ... U_#cutoff)
				Matrix cov_U,cov_U_tmp;
				if(condition)
					cov_U_tmp.Product(cov_weight,cond_cov[g]);
				else
					cov_U_tmp.Product(cov_weight,cov[g]);
				Matrix cov_weight_trans;
				cov_weight_trans.Transpose(cov_weight);
				cov_U.Product(cov_U_tmp,cov_weight_trans); //now, cov(U) is saved in cov_U
				//Calculate covariance matrix for (T_1 ... T_#cutoff)
				Matrix cov_T;
				cov_T.Dimension(cov_U.rows,cov_U.cols);
				cov2cor(cov_U,cov_T);
				
				pvalue = CalculateMVTPvalue(score,cov_T,t_max);
				printf("completed!\n");
			}
			else {
				pvalue = CalculateMVTPvalue(score,cov_T,t_max);
				printf("completed!\n");
			}
		}
		else
			pvalue = 1.0 - result[0];
		
		if((condition && cond!="") || cond=="")
		{
			if(pvalue <report_pvalue_cutoff && report)
			{
				StringArray variants;
				variants.AddTokens(var,";");
				for(int v=0;v<tmp_maf.Length();v++)
					ifprintf(reportOutput,"%s\t%s\t%g\t%g\t%g\t%s\t%g\t%g\t%g\n",group.annoGroups[g].c_str(),method.c_str(),pvalue,MAF_cutoff,chosen_cutoff,variants[v].c_str(),tmp_maf[v],tmp_eff[v],tmp_pvalue[v]);
			}
		}
		
		if(cond=="" || (!condition && cond!=""))
		{
			if(fullResult)
			{
				ifprintf(output,"%s\t%d\t%s\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str());
				
				for(int i=0;i<tmp_maf.Length()-1;i++)
					ifprintf(output,"%g,",tmp_maf[i]);
				ifprintf(output,"%g\t",tmp_maf[tmp_maf.Length()-1]);
				
				for(int i=0;i<tmp_eff.Length()-1;i++)
					ifprintf(output,"%g,",tmp_eff[i]);
				ifprintf(output,"%g\t",tmp_eff[tmp_eff.Length()-1]);
				
				for(int i=0;i<tmp_pvalue.Length()-1;i++)
					ifprintf(output,"%g,",tmp_pvalue[i]);
				ifprintf(output,"%g\t",tmp_pvalue[tmp_pvalue.Length()-1]);
				
				ifprintf(output,"%g\t%g\t%g\t%g\t%g\t%g\t",average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
			}
			else
				ifprintf(output,"%s\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t",group.annoGroups[g].c_str(),tmp_maf.Length(),var.c_str(),average_af,min_af,max_af,chosen_effSize,chosen_cutoff,pvalue);
			if(cond=="")
				ifprintf(output,"\n");
		}
		
		if(cond!="" && condition)
			ifprintf(output,"%g\t%g\t%g\n",chosen_effSize,chosen_cutoff,pvalue);
		
		if(pvalue>1.0)
		pvalue = 1.0;
	}
	if(lower) delete [] lower;
	if(upper) delete [] upper;
	if(mean) delete [] mean;
	return pvalue;
}


void annoGroups::printVTheader(IFILE& output)
{
	ifprintf(output,"##Method=VT\n");
	ifprintf(output,"##STUDY_NUM=%d\n",scorefile.Length());
	ifprintf(output,"##TotalSampleSize=%d\n",total_N);
	if(fullResult)
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tMAFs\tSINGLEVAR_EFFECTs\tSINGLEVAR_PVALUEs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
	else
		ifprintf(output,"#GROUPNAME\tNUM_VAR\tVARs\tAVG_AF\tMIN_AF\tMAX_AF\tEFFECT_SIZE\tMAF_CUTOFF\tPVALUE\t");
	
	if(cond!="")
		ifprintf(output,"EFFECT_SIZE\tMAF_CUTOFF\tCOND_PVALUE\n");
	else
		ifprintf(output,"\n");
}
*/