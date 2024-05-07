package gov.fda.gsrs.ndsri;

public class FeatureScorePair{
	
	String featureName;
	
	String scoreName;
	
	public FeatureScorePair(String featureName,String scoreName){
		this.featureName=featureName;
		this.scoreName=scoreName;			
	}
	
	public FeatureScorePairInstance getInstance(String value, int score){
		return new FeatureScorePairInstance(this, value,score);
	}
	
	public FeatureScorePairInstance getInstanceYesNo(String value, int scoreIfYes){
		if(value.equals("YES")){
			return new FeatureScorePairInstance(this, value,scoreIfYes);
		}else{
			return new FeatureScorePairInstance(this, value,0);
		}
	}
	public String getFeatureName() {
		return featureName;
	}

	public void setFeatureName(String featureName) {
		this.featureName = featureName;
	}

	public String getScoreName() {
		return scoreName;
	}

	public void setScoreName(String scoreName) {
		this.scoreName = scoreName;
	}

}