package uk.ac.cam.ch.wwmm.opsin;
/**
 * 
 * @author dl387
 *
 * @param <T>
 */
class PropertyKey<T> {
    private final String name;

    public PropertyKey(String name) {
        this.name = name;
    }

    @Override
    public int hashCode() {
        return 37 * (name != null ? name.hashCode() : 0);
    }

    @Override
    public boolean equals(Object obj) {
        return this == obj;
    }

    @Override
    public String toString() {
        return "Key{" + name + "}";
    }
}
