package gov.fda.gsrs.chem.util;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.fda.gsrs.chem.util.StreamUtil.Final;
import gov.nih.ncats.common.Tuple;
import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;


public class AtomTest extends ChemWrapper<Atom>{
	private final Atom ca;
	
	public static Stream<AtomTest> stream(Chemical c){
		
		
		return c.atoms()
				     .map(AtomTest::of);
	}
	
	public AtomTest(Atom ca2) {
		this.ca=ca2;
	}
	public boolean hasAromaticBond() {
		
		return ca.hasAromaticBond();
	}

	public String getAtomSymbol(){
		return ca.getSymbol();
	}
	
	public boolean isAtomSymbol(String s){
		return ca.getSymbol().equals(s);
	}
	
	public boolean isRingAtom(){
		return ca.isInRing();
	}
	public boolean isInRing(){
		return ca.isInRing();
		
	}
	
	public boolean isQueryAtom(){
		return ca.isQueryAtom();
	}
	
	/**
	 * Is N of the form RN(C=O)R'
	 * @return
	 */
	public boolean isPeptideAmine(){
		return isAtomSymbol("N") && hasNeighbor((cn)->cn.isCarbonyl() && cn.hasCarbon()) 
				                 && hasNeighbor((cn)->cn.isCarbon() && !cn.isCarbonyl());
	}
	
	
	
	
	
	public boolean isNitrosamineAtom(){
		return     this.isNitrogen()
				&& this.isNeutral()
				&& this.getBonds().count()==2
				&& this.hasNeighbor(at->at.isOxygen() && at.hasDoubleBond())
				&& this.hasNeighbor(at->at.isNitrogen());
	}
	
	public boolean hasCarbon(){
		return hasNeighbor((ct)->ct.isCarbon());
	}
	public int countCarbon(){
		return countNeighbor((ct)->ct.isCarbon());
	}

	public Stream<BondTest> getBonds(){
		return ca.getBonds().stream().map(BondTest::of);
	}
	
	public Stream<AtomTest> getNeighbors(){
		return ca.getNeighbors()
				.stream()
				.map(AtomTest::of);
	}
	
	public LinkedList<AtomTest> getFirstPathTo(Predicate<AtomTest> at){
	
		Final<LinkedList<AtomTest>> fll = Final.of(null);
		
		getBreadthFirstPathsUntil(l->{
			AtomTest atom=l.peek();
			if(at.test(atom)){
				fll.set(l);
				return false;
			}
			return true;
		});
		return fll.get();
	}
	
	public Optional<BondTest> getBondTo(AtomTest oat){
		return this.getBonds().filter(b->b.hasAtom(oat)).findFirst();
	}
	
	public LinkedList<AtomTest> getFirstPathTo(AtomTest otherAtom){
		return getFirstPathTo(cn->otherAtom.equals(cn));
	}
		
	
//	public AtomTest getFirstChain
	
	public AtomTest getBreadthFirstPathsUntil(Function<LinkedList<AtomTest>,Boolean> cont, Set<AtomTest> excludeThru){
		Set<AtomTest> already = new HashSet<>();
		already.add(this);
		already.addAll(excludeThru);
		
		
		Set<AtomTest> currentLayer = new HashSet<>();
		currentLayer.add(this);
		
		LinkedList<AtomTest> path = new LinkedList<>();
		path.push(this);
		if(!cont.apply(path))return this;
		
		Map<AtomTest, LinkedList<AtomTest>> as= new HashMap<>();
		
		as.put(this, path);
		
		
		while(!currentLayer.isEmpty()){
			
			Set<AtomTest> ncurrentLayer =  currentLayer.stream()
					    .map(ca->Tuple.of(ca,as.get(ca)))
			            .flatMap(t->t.k().getNeighbors()
			            		         .peek(ca->{
			            		        	 LinkedList<AtomTest> allist=as.computeIfAbsent(ca, (k)->new LinkedList<>(t.v()));
			            		        	 if(!allist.contains(ca)){
			            		        		 allist.push(ca);
			            		        	 }
			            		         })
			            				)
			            .filter(ca->!already.contains(ca))
			            .collect(Collectors.toSet());
			
			for(AtomTest at:ncurrentLayer){
				if(!cont.apply(as.get(at))){
					return this;
				}
			}
			
			already.addAll(ncurrentLayer);
			currentLayer=ncurrentLayer;
			
		}
		return this;
	}
	
	public AtomTest getBreadthFirstPathsUntil(Function<LinkedList<AtomTest>,Boolean> cont){
		Set<AtomTest> empty= new HashSet<>();
		
		return getBreadthFirstPathsUntil(cont, empty);
	}
		
		
	
	public boolean isCarbon(){
		return isAtomSymbol("C");
	}
	
	public boolean isHalogen(){
		return isAtomSymbol("Cl") || isAtomSymbol("F") || isAtomSymbol("Br") || isAtomSymbol("I");
	}
	
	
	public boolean isHydrogen(){
		return isAtomSymbol("H");
	}
	
	public boolean isOxygen(){
		return isAtomSymbol("O");
	}
	public boolean isNitrogen(){
		return isAtomSymbol("N");
	}
	
	
	
	public int getDistanceTo(AtomTest ca2){
		
		if(this.equals(ca2)){
			return 0;
		}
		int ndist = 1;
		
		Set<AtomTest> alreadySaw = new HashSet<AtomTest>();
		alreadySaw.add(this);
		Set<AtomTest> cNeighbors = new HashSet<AtomTest>();
		cNeighbors.addAll(this.getNeighbors().collect(Collectors.toList()));
		while(ndist<10){
			if(cNeighbors.contains(ca2)){
				return ndist;
			}else{
				alreadySaw.addAll(cNeighbors);			
				cNeighbors = cNeighbors.stream()
				          .flatMap(caa->caa.getNeighbors())
				          .filter(caa->!alreadySaw.contains(caa))
				          .collect(Collectors.toSet());
				ndist++;
			}
		}
		return Integer.MAX_VALUE;
	}
	
	
	
	
	public boolean isCarbonChainOH(){
		return isOH() && hasNeighbor(ca2->ca2.isCarbon() && !ca2.isInRing());
	}
	
	public boolean isCarbonInRingConnectedToNitrogenInRingAndOxygenInRing(){
		return isCarbon() && isInRing() 
				&& hasNeighbor(ca2->ca2.isNitrogen() && ca2.isRingAtom())
				&& hasNeighbor(ca2->ca2.isOxygen() && ca2.isRingAtom());
	}
	
	public boolean isLikelyFivePrime(){
		return isCarbonChainOH()
				&& this.hasNeighbor(ca2->ca2.hasNeighbor(ca3->ca3.isLikelyFivePrimeRingCarbon()));
	}
	
	public boolean isLikelyThreePrime(){
		return this.isCarbonRingOH() 
				&& !this.hasNeighbor(ca2->ca2.hasNeighbor(ca3->ca3.isCarbonInRingConnectedToNitrogenInRingAndOxygenInRing()))
				&& this.hasNeighbor(ca2->ca2.hasNeighbor(ca3->ca3.isLikelyTwoPrimeCarbon()));
	}
	
	public boolean isLikelyTwoPrimeCarbon(){
		return this.isCarbon() && this.isRingAtom() &&
			   this.hasNeighbor(ca2->ca2.isLikelyNucleoBaseSite());
	}
	
	public boolean isLikelyOxygenInSugar(){
		return this.isOxygenInRing() && this.hasNeighbor(ca2->ca2.isLikelyNucleoBaseSite());
	}
	
	public boolean isLikelyFivePrimeRingCarbon(){
		return this.isCarbon() && this.isRingAtom() && this.hasNeighbor(ca2->ca2.isLikelyOxygenInSugar())
			   && !this.isLikelyNucleoBaseSite();
	}
	
	
	
	public boolean isCarbonRingOH(){
		return isOH() && hasNeighbor(ca2->ca2.isCarbon() && ca2.isInRing());
	}
	
	public boolean isOxygenInRing(){
		return isOxygen() && getHCount()==0 && isInRing();
	}
	
	public boolean isCarbonRingConnectedNitrogen(){
		
		return isNitrogen() && attachedToCarbonRing();
	}
	
	public boolean attachedToCarbonRing(){
		return hasNeighbor(ca2->ca2.isCarbon() && ca2.isInRing());
	}
	
	public boolean isLikelyNucleoBaseSite(){
		return this.isCarbon() && this.hasNeighbor(ca2->ca2.isOxygenInRing()) && this.hasNeighbor(ca2->ca2.isCarbonRingConnectedNitrogen());
	}
	
	public boolean isPhosphorus(){
		return isAtomSymbol("P");
	}
	public boolean isSulfur(){
		return isAtomSymbol("S");
	}
	
	public Stream<AtomTest> getNeighbors(Predicate<AtomTest> atomtest){
		return getNeighbors().filter(atomtest);
	}
	
	public boolean hasNeighbor(Predicate<AtomTest> atomtest){
		return getNeighbors(atomtest).findAny().isPresent();
	}
	
	
	public boolean hasBond(Predicate<BondTest> bt){
		return getBonds().anyMatch(bt);
	}
	
	public boolean hasDoubleBond(){
		return hasBond((bt)->bt.isDoubleBond());
	}
	public boolean hasDoubleBondCarbon(){
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isDoubleBondedCarbon() && cn.hasBond(b->b.hasAtom(this) && b.isDoubleBond());
		});
	}
	
	public boolean hasCarbonTripleBond(){
		return isCarbon() && this.hasBond(bt->bt.isTripleBond());
	}
	
	public boolean isDoubleBondedOxygen(){
		return isOxygen() && hasDoubleBond();
	}

	public boolean isDoubleBondedCarbon(){
		return isCarbon() && hasDoubleBond();
	}
	
	public boolean isDoubleBondedNitrogen(){
		return isNitrogen() && hasDoubleBond();
	}
	
	public boolean isDoubleBondedNSorP(){
		return (isSulfur() || isNitrogen() || isPhosphorus()) && (hasDoubleBond());
	}
	
	public boolean isDoubleBondedOrAromaticNSorP(){
		
		if(isNitrogen() && hasAromaticBond()) {
			return getHCount()==0;
		}
		
		return (isSulfur() || isNitrogen() || isPhosphorus()) && (hasDoubleBond() || hasAromaticBond())
				;
	}
	
	public boolean isCarbonyl(){
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isDoubleBondedOxygen();
		});
	}
	
	public boolean isCarbonylAmide(){
		return isCarbonyl() && this.hasBond(b->b.isSingleBond() && b.hasAtomKind(att->att.isNitrogen()));
	}
	
	public boolean isCarbonylEster(){
		return isCarbonyl() && this.hasBond(b->b.isSingleBond() && b.hasAtomKind(att->att.isOxygen() && !att.isOH()));
	}

	public boolean isImine(){
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isDoubleBondedNitrogen() && cn.hasBond(b->b.hasAtom(this) && b.isDoubleBond());
		});
	}
	public boolean isSulfone(){
		return isSulfur() && countNeighbor((cn)->{
			return cn.isDoubleBondedOxygen();
		})==2;
	}
	
	public boolean isNitrate(){
		return isNitrogen() && this.getHCount()==0 && countNeighbor((cn)->{
			return cn.isDoubleBondedOxygen();
		})==1 && this.hasOH();
	}
	
	public boolean isQuatAmine(){
		return isNitrogen() && this.getHCount()==0 && this.getBonds().count()==4;
	}
	
	
	public boolean isPhosphone(){
		return isPhosphorus() && countNeighbor((cn)->{
			return cn.isDoubleBondedOxygen();
		}) >=1;
	}
	
	public boolean isEnolLike(){
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isOH();
		}) && (this.hasDoubleBond() || this.isAromaticCarbon());
	}
	
	public boolean isTautomerOfDoubleBondedNSorPIgnoring(AtomTest at2){
		return isCarbon() && hasNeighbor((cn)->{
			if(at2.equals(cn))return false;
			return (cn.isSulfur() || cn.isNitrogen() || cn.isPhosphorus()) && cn.getBondTo(this).map(b->b.isSingleBond()).isPresent();
		}) && (this.hasDoubleBond() || this.isAromaticCarbon());
	}
	
	public boolean isMethyl(){
		return isCarbon() && this.getHCount()==3;
	}
	
	public boolean isTerminalHeteroAtom(){
		return !isCarbon() && (this.countNeighbor(aa->!aa.isHydrogen())<=1);
	}
	
	public boolean isCarboxyl(){
		return isCarbonyl() && hasOH();
	}
	public boolean isCarboxylCharged(){
		return isCarbonyl() && hasNeighbor((cn)->cn.isOxygen() && cn.getHCount()==0 && !cn.hasDoubleBond() && !cn.isNeutral() );
	}
	
	public boolean isOH(){
		return isOxygen() && getHCount()>0;
	}
	public boolean hasOH(){
		return hasNeighbor((cn)->cn.isOH());
	}
	
	public int countOH(){
		return countNeighbor((cn)->cn.isOH());
	}
	
	public int getHCount(){
		return ca.getImplicitHCount() + (int)getNeighbors(cn->cn.isAtomSymbol("H")).count();
	}
	
	public boolean isSp3Carbon(){
		
		return isCarbon() && (ca.getImplicitHCount() + ca.getBondCount()==4);
	}
	public boolean hasAllCarbonAndHNeighbors(){
		return !hasNeighbor(nn->!nn.isCarbon() && !nn.get().getSymbol().equals("H"));
	}
	public boolean isAromaticCarbon(){
		return isCarbon() && (ca.hasAromaticBond()) && isInRing();
	}
	
	public boolean isAromatic(){
		return ca.hasAromaticBond() && isInRing();
	}
	
	@Override
	public Atom get() {
		return ca;
	}

	public static AtomTest of(Atom ca){
		return new AtomTest(ca);
	}
	
	
	public boolean equals(Object o){
		if(!(o instanceof AtomTest)){
			return false;
		}
		return this.wEquals((AtomTest)o);
	}
	public int hashCode(){
		return this.wHashCode();
	}

	public int countNeighbor(Predicate<AtomTest> atomtest) {
		return (int)getNeighbors(atomtest).count();
	}

	public boolean isDoubleBondToSNorP() {
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isDoubleBondedNSorP() && cn.hasBond(b->b.hasAtom(this) && b.isDoubleBond());
		});
	}
	
	public boolean isDoubleOrAromaticBondToSNorP() {
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isDoubleBondedOrAromaticNSorP() && cn.hasBond(b->b.hasAtom(this) && (b.isDoubleBond()|| b.isAromatic()));
		});
	}
	
	public boolean isDoubleBondToSorN() {
		return isCarbon() && hasNeighbor((cn)->{
			return !cn.isPhosphorus() && cn.isDoubleBondedNSorP() && cn.hasBond(b->b.hasAtom(this) && b.isDoubleBond());
		});
	}
	
	public boolean isCyanide() {
		return isCarbon() && hasNeighbor((cn)->{
			return cn.isNitrogen() && cn.hasBond(b->b.hasAtom(this) && b.isTripleBond());
		});
	}

	public boolean isNeutral() {
		return get().getCharge()==0;
	}

	public boolean isSulfoxyl(){
		return isSulfone() && hasOH();
	}
	
	public boolean isSulfoxylCharged(){
		return isSulfone() && hasNeighbor((cn)->cn.isOxygen() && cn.getHCount()==0 && !cn.hasDoubleBond()  && !cn.isNeutral());
	}
	

	public boolean  isPhosphate() {
		return isPhosphone() && countOH()>=2;
	}

	public boolean isMetal() {
		return ca.isMetal();
	}
	
}