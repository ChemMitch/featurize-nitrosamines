package gov.fda.gsrs.chem.util;


public abstract class ChemWrapper<T> implements Wrapper<T>{
	public boolean equals(Object o){
		if(o==null)return false;
		if(!(o instanceof Wrapper))return false;
		return this.wEquals((Wrapper<T>)o);
	}
	public int hashCode(){
		return this.wHashCode();
	}
}