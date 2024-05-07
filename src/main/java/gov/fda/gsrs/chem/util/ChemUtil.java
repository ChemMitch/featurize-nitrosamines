package gov.fda.gsrs.chem.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;

public class ChemUtil {
	
	public static void simpleCleanup(Chemical c){
		//precleaning needed for bad metal reading
		AtomTest.stream(c)
		.filter(at->at.isMetal())
		.forEach(at->{
			at.get().setImplicitHCount(0);
		});
	}
	
	public static void addChemical(Chemical parent, Chemical child) {
		Map<Atom,Atom> newAtomMap = new HashMap<>();
		
		child.atoms().forEach(atOld->{
			
			Atom atNew = parent.addAtom(atOld.getSymbol());
			atOld.getAlias().ifPresent(ali->{
				atNew.setAlias(ali);
			});
			atOld.getAtomToAtomMap().ifPresent(ali->{
				atNew.setAtomToAtomMap(ali);
			});
			atNew.setAtomCoordinates(atOld.getAtomCoordinates());
			atNew.setCharge(atOld.getCharge());
			atNew.setMassNumber(atOld.getMassNumber());
			
			atNew.setRadical(atOld.getRadical());
			if(atOld.getRGroupIndex().isPresent()) {
				atNew.setRGroup(atOld.getRGroupIndex().getAsInt());				
			}
			newAtomMap.put(atOld, atNew);
		});
		child.bonds().forEach(bOld->{
			Atom nat1 = newAtomMap.get(bOld.getAtom1());
			Atom nat2 = newAtomMap.get(bOld.getAtom2());
			Bond bNew=parent.addBond(nat1,nat2,bOld.getBondType());
			bNew.setStereo(bOld.getStereo());
		});
		child.atoms().forEach(atOld->{
			Atom atNew = newAtomMap.get(atOld);
			atNew.setImplicitHCount(atOld.getImplicitHCount());
		});
	}

	
	public static List<Bond> getNeighborBonds(Chemical c, Bond cb){
		List<Bond> bonds = new ArrayList<>();
		
		for(Atom at:new Atom[] {cb.getAtom1(),cb.getAtom2()} ){
			at.getBonds()
			  .forEach(cbb->{
				  if(!cbb.equals(cb)){
					  bonds.add(cbb);
				  }
			  });
		}
		return bonds;
	}
	
	/**
	 * This is a ring-detection algorithm which colors edges based on their smallest detected ring, up to the maxRingSize
	 * specified. This is accomplished with a breadth-first search from each bond
	 * @param b
	 * @param m
	 * @return
	 */
	public static Map<Integer,Integer> getSmallestRingSizeForEachBond(Chemical b, int maxRingSize){

		//The computed return object
		Map<Integer,Integer> minRingForBond= new HashMap<>();
		
		// The neighbor bonds "adjacency" type map. Each bond points
		// to the set of bonds it shares an atom with
		Map<Integer,BitSet> nbonds= new HashMap<>();
		
		
		// Lookups to go from a bond or atom to its index
		// with some effort these could be eliminated
		Map<Bond, Integer> bondIndex = new HashMap<>();
		Map<Atom, Integer> atomIndex = new HashMap<>();
		
		//Cursor Bonds (left and right)
		Map<Integer,Set<Integer>> cursorBondsL= new HashMap<>();
		Map<Integer,Set<Integer>> cursorBondsR= new HashMap<>();
		
		//Cursor Atoms (left and right)
		Map<Integer,Set<Integer>> cursorAtomsL= new HashMap<>();
		Map<Integer,Set<Integer>> cursorAtomsR= new HashMap<>();
		Map<Integer,Set<Integer>> previousBondsL= new HashMap<>();
		Map<Integer,Set<Integer>> previousBondsR= new HashMap<>();
		Map<Integer,Set<Integer>> previousAtomsL= new HashMap<>();
		Map<Integer,Set<Integer>> previousAtomsR= new HashMap<>();
		
		
		for(int i=0;i<b.getAtomCount();i++){
			atomIndex.put(b.getAtom(i), i);
		}
		for(int i=0;i<b.getBondCount();i++){
			Bond cb=b.getBond(i);
			bondIndex.put(cb, i);
		}
		for(int i=0;i<b.getBondCount();i++){
			Bond cb=b.getBond(i);
//			if(cb.isRingBond()){
				previousBondsL.put(i,  Stream.of(i).collect(Collectors.toSet()));
				previousBondsR.put(i,  Stream.of(i).collect(Collectors.toSet()));
				
				previousAtomsL.put(i,  Stream.of(atomIndex.get(cb.getAtom1())).collect(Collectors.toSet()));
				previousAtomsR.put(i,  Stream.of(atomIndex.get(cb.getAtom2())).collect(Collectors.toSet()));
				
				cursorBondsL.put(i, cb.getAtom1().getBonds().stream().filter(bbb->!bbb.equals(cb)).map(cbb->bondIndex.get(cbb)).collect(Collectors.toSet()));
				cursorBondsR.put(i, cb.getAtom2().getBonds().stream().filter(bbb->!bbb.equals(cb)).map(cbb->bondIndex.get(cbb)).collect(Collectors.toSet()));
				
				Set<Integer> patsl= previousAtomsL.get(i);
				Set<Integer> patsr= previousAtomsR.get(i);
				
				cursorAtomsL.put(i, cursorBondsL.get(i).stream()
						                               .flatMap(bbi->Arrays.stream(getAtoms(b.getBond(bbi))))
						                               .map(at->atomIndex.get(at))
						                               .filter(ati->!patsl.contains(ati))
						                               .collect(Collectors.toSet())
						                               );
				cursorAtomsR.put(i, cursorBondsR.get(i).stream()
                        .flatMap(bbi->Arrays.stream(getAtoms(b.getBond(bbi))))
                        .map(at->atomIndex.get(at))
                        .filter(ati->!patsr.contains(ati))
                        .collect(Collectors.toSet())
                        );
//			}
		}
		for(int i=0;i<b.getBondCount();i++){
			Bond cb=b.getBond(i);
			BitSet bsBond=new BitSet(b.getBondCount());
			getNeighborBonds(b, cb)
				.stream()
				.map(cbb-> bondIndex.get(cbb))
				.filter(Objects::nonNull)
				.forEach(ii->{
					bsBond.set(ii);
				});
			nbonds.put(i,bsBond);
		}
		int[] cii= new int[]{0,0};
		int maxIter = (int) Math.ceil((maxRingSize-3)*0.5);
		for(int iter=0;iter<=maxIter;iter++){
			cii[1]=iter;
			for(int i=0;i<b.getBondCount();i++){
				cii[0]=i;

				Set<Integer> cbl= cursorBondsL.get(i);
				Set<Integer> cbr= cursorBondsR.get(i);
				Set<Integer> casl= cursorAtomsL.get(i);
				Set<Integer> casr= cursorAtomsR.get(i);
				
				if(minRingForBond.get(i)==null && cbl!=null && cbr!=null){
					boolean ringBonds= cbl.stream().filter(cali->cbr.contains(cali)).findAny().isPresent();
					if(ringBonds){
						minRingForBond.put(i, iter*2+2);
						continue;
					}
					if(iter*2+3<=maxRingSize){
						boolean ringatoms= casl.stream().filter(cali->casr.contains(cali)).findAny().isPresent();
						if(ringatoms){
							minRingForBond.put(i, iter*2+3); //iter*2+3 = maxRing  iter=(maxRing-3)/2
							continue;
						}
					}
					
					Set<Integer> pbr= previousBondsR.computeIfAbsent(i, k->{
						return new HashSet<>();
					});
					Set<Integer> pbl= previousBondsL.computeIfAbsent(i, k->{
						return new HashSet<>();
					});
					Set<Integer> par= previousAtomsR.computeIfAbsent(i, k->{
						return new HashSet<>();
					});
					Set<Integer> pal= previousAtomsL.computeIfAbsent(i, k->{
						return new HashSet<>();
					});
					Set<Integer> newCursorR = new HashSet<>();
					Set<Integer> newCursorL = new HashSet<>();
					Set<Integer> newAtomsR  = new HashSet<>();
					Set<Integer> newAtomsL  = new HashSet<>();
					
					cbl.stream()
					  .flatMap(j->nbonds.get(j).stream().mapToObj(ib->(Integer)ib))
					  .filter(bi->!pbl.contains(bi))
					  .filter(bi->!cbl.contains(bi))
					  .forEach(ii->{
						  newCursorL.add(ii);
					  });
					cbr.stream()
					  .flatMap(j->nbonds.get(j).stream().mapToObj(ib->(Integer)ib))
					  .filter(bi->!pbr.contains(bi))
					  .filter(bi->!cbl.contains(bi))
					  .forEach(ii->{
						  newCursorR.add(ii);
					  });
					
					newCursorL.stream()
					         .flatMap(ii->Arrays.stream(getAtoms(b.getBond(ii))))
					         .map(ca->atomIndex.get(ca))
					         .filter(ii->!casl.contains(ii))
					         .forEach(ii->{
								  newAtomsL.add(ii);
							  });
					newCursorR.stream()
					          .flatMap(ii->Arrays.stream(getAtoms(b.getBond(ii))))
					          .map(ca->atomIndex.get(ca))
					          .filter(ii->!casr.contains(ii))
					          .forEach(ii->{
						 		  newAtomsR.add(ii);
					  		  });
					         
					pbr.addAll(cbr);
					pbl.addAll(cbl);
					par.addAll(casr);
					pal.addAll(casl);
					
					cursorBondsR.put(i,newCursorR);
					cursorBondsL.put(i,newCursorL);
					cursorAtomsR.put(i,newAtomsR);
					cursorAtomsL.put(i,newAtomsL);
					
				}
			}
		}
		
		return minRingForBond;
	}
	
	public static Atom[] getAtoms(Bond b){
		return new Atom[] {b.getAtom1(),b.getAtom2()};
	}
	
	/**
	 * BitSet is 1-indexed. Returns cannonical key for chemical, with
	 * given atoms at specified indices removed. This is useful
	 * for determining symmetry.
	 * 
	 * @param cc
	 * @param atoms
	 * @return
	 */
	public static String getHashWithout(Chemical cc, MuteBitSet atoms){
		Chemical c= cc.copy();
		
		Set<Atom> toRemove = atoms.stream()
		     .filter(i->i>0 && i<=c.getAtomCount())
		     .mapToObj(i->c.getAtom(i-1))
		     .collect(Collectors.toSet());
		
		toRemove.forEach(ca->{
			c.removeAtom(ca);
		});
		
		try {
			return c.toInchi().getKey();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		return null;
		
	}
}
