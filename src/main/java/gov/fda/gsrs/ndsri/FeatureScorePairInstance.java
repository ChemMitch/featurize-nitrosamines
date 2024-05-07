package gov.fda.gsrs.ndsri;

public class FeatureScorePairInstance{
	FeatureScorePair pair;
	

	String value;
	int score;
	
	public FeatureScorePairInstance(FeatureScorePair pair, String value, int score){
		this.pair=pair;
		this.value=value;
		this.score=score;
		
	}
	
	public FeatureScorePair getPair() {
		return pair;
	}

	public void setPair(FeatureScorePair pair) {
		this.pair = pair;
	}

	public String getValue() {
		return value;
	}

	public void setValue(String value) {
		this.value = value;
	}

	public int getScore() {
		return score;
	}

	public void setScore(int score) {
		this.score = score;
	}
	
}