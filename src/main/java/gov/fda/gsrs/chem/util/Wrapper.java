package gov.fda.gsrs.chem.util;

public interface Wrapper<T>{
	public T get();
	
	public default int wHashCode(){
		T t=get();
		if(t==null)return 0;
		return t.hashCode();
	}
	
	public default boolean wEquals(Wrapper<T> w){
		T t=get();
		T t2=w.get();
		if(t==null)return (t2==null);
		return t.equals(t2);
	}
}