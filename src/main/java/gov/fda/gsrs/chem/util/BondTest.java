package gov.fda.gsrs.chem.util;

import java.util.function.Predicate;
import java.util.stream.Stream;

import gov.nih.ncats.common.stream.StreamUtil;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;


public class BondTest extends ChemWrapper<Bond>{
	private final Bond cb;
	
	public static Stream<BondTest> stream(Chemical c){
		return c.bonds()
				     .map(BondTest::of);
	}
	
	public BondTest(Bond cb){
		this.cb=cb;
	}
	
	public boolean isInRing() {
		return cb.isInRing();
	}
	
	public boolean isType(int bondType){
		return cb.getBondType().ordinal() == bondType;
	}

	public boolean hasAtom(AtomTest at){
		return cb.getAtom1().equals(at.get()) || cb.getAtom2().equals(at.get());
	}
	
	public boolean hasAtomKind(Predicate<AtomTest> pred){
		return pred.test(AtomTest.of(cb.getAtom1())) || pred.test(AtomTest.of(cb.getAtom2()));
	}
	
	@Override
	public Bond get() {
		return cb;
	}
	
	public boolean isDoubleBond(){
		return cb.getBondType().getOrder()==2;
	}
	
	public boolean isSingleBond(){
		return cb.getBondType().getOrder()==1;
	}
	public boolean isTripleBond(){
		return cb.getBondType().getOrder()==3;
	}
	
	public boolean isAromatic() {
		return cb.isAromatic();
	}
	public Stream<AtomTest> atoms(){
		return Stream.of(this.get().getAtom1(), this.get().getAtom2()).map(a->AtomTest.of(a));
	}
	
	public Stream<BondTest> getNeighborBonds(){
		return StreamUtil.with(this.get().getAtom1().getBonds().stream().map(b->BondTest.of(b)))
				.and(this.get().getAtom2().getBonds().stream().map(b->BondTest.of(b)))
				.stream()
				.distinct()
				.filter(bt->!this.equals(bt))
				;
	}
	
	public static BondTest of(Bond cb){
		return new BondTest(cb);
	}
	
	public String getType(){
		if(isSingleBond())return"-";
		if(isDoubleBond())return"=";
		if(isTripleBond())return"#";
		return"~";
	}
	

	public boolean equals(Object o){
		if(!(o instanceof BondTest)){
			return false;
		}
		return this.wEquals((BondTest)o);
	}
	public int hashCode(){
		return this.wHashCode();
	}
}