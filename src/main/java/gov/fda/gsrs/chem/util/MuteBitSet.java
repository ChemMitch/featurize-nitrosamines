package gov.fda.gsrs.chem.util;

import java.util.BitSet;
import java.util.Set;
import java.util.stream.IntStream;

public class MuteBitSet{
	private BitSet internal= new BitSet();
	
	public MuteBitSet(){
		
	}
	
	public MuteBitSet(BitSet bs){
		this.internal=bs;
	}

	public byte[] toByteArray() {
		return internal.toByteArray();
	}

	public long[] toLongArray() {
		return internal.toLongArray();
	}

	public MuteBitSet flip(int bitIndex) {
		internal.flip(bitIndex);
		return this;
	}

	public MuteBitSet flip(int fromIndex, int toIndex) {
		internal.flip(fromIndex, toIndex);
		return this;
	}

	public MuteBitSet set(int bitIndex) {
		internal.set(bitIndex);
		return this;
	}

	public MuteBitSet set(int bitIndex, boolean value) {
		internal.set(bitIndex, value);
		return this;
	}

	public MuteBitSet set(int fromIndex, int toIndex) {
		internal.set(fromIndex, toIndex);
		return this;
	}

	public MuteBitSet set(int fromIndex, int toIndex, boolean value) {
		internal.set(fromIndex, toIndex, value);
		return this;
	}

	public MuteBitSet clear(int bitIndex) {
		internal.clear(bitIndex);
		return this;
	}

	public MuteBitSet clear(int fromIndex, int toIndex) {
		internal.clear(fromIndex, toIndex);
		return this;
	}

	public MuteBitSet clear() {
		internal.clear();
		return this;
	}

	public boolean get(int bitIndex) {
		return internal.get(bitIndex);
	}

	public MuteBitSet get(int fromIndex, int toIndex) {
		return new MuteBitSet(internal.get(fromIndex, toIndex));
	}

	public int nextSetBit(int fromIndex) {
		return internal.nextSetBit(fromIndex);
	}

	public int nextClearBit(int fromIndex) {
		return internal.nextClearBit(fromIndex);
	}

	public int previousSetBit(int fromIndex) {
		return internal.previousSetBit(fromIndex);
	}

	public int previousClearBit(int fromIndex) {
		return internal.previousClearBit(fromIndex);
	}

	public int length() {
		return internal.length();
	}

	public boolean isEmpty() {
		return internal.isEmpty();
	}

	public boolean intersects(BitSet set) {
		return internal.intersects(set);
	}

	public int cardinality() {
		return internal.cardinality();
	}

	public MuteBitSet and(MuteBitSet muteBitSet) {
		internal.and(muteBitSet.toBitSet());
		return this;
	}

	public MuteBitSet or(MuteBitSet set) {
		internal.or(set.toBitSet());
		return this;
	}

	public MuteBitSet xor(MuteBitSet set) {
		internal.xor(set.toBitSet());
		return this;
	}

	public MuteBitSet andNot(MuteBitSet set) {
		internal.andNot(set.toBitSet());
		return this;
	}
	


	
	
	
	public MuteBitSet newAnd(MuteBitSet muteBitSet) {
		return new MuteBitSet((BitSet)internal.clone()).and(muteBitSet);
	}

	public MuteBitSet newOr(MuteBitSet set) {
		return new MuteBitSet((BitSet)internal.clone()).or(set);
	}

	public MuteBitSet newXor(MuteBitSet set) {
		return new MuteBitSet((BitSet)internal.clone()).xor(set);
	}
	
	public MuteBitSet newAndNot(MuteBitSet set) {
		return new MuteBitSet((BitSet)internal.clone()).andNot(set);
	}
	
	
	

	public int hashCode() {
		return internal.hashCode();
	}

	public int size() {
		return internal.size();
	}

	public boolean equals(Object obj) {
		return internal.equals(obj);
	}

	public Object clone() {
		return internal.clone();
	}

	public String toString() {
		return internal.toString();
	}

	public IntStream stream() {
		return internal.stream();
	}
	
	public BitSet toBitSet(){
		return this.internal;
	}
	
	
	public static MuteBitSet from(Set<Integer> iset){
		MuteBitSet bs = new MuteBitSet();
		iset.stream().forEach(i->bs.set(i));
		return bs;
	}
	
	public static MuteBitSet from(int ...is){
		MuteBitSet bs = new MuteBitSet();
		IntStream.of(is)
		.forEach(i->{
			bs.set(i);
		});
		return bs;
	}
}