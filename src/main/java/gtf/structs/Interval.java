package gtf.structs;

import java.util.Objects;

public class Interval implements Comparable<Interval> {
    // one-based, inclusive
    private int start;
    private int end;

    public Interval(int start, int end) {
        this.start = start;
        this.end = end;

    }

    public Interval(Interval interval) {
        this.start = interval.start;
        this.end = interval.end;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Interval interval = (Interval) o;
        return start == interval.start && end == interval.end;
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end);
    }

    @Override
    public String toString() {
        return start + ":" + (end + 1);
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getLength() {
        return end - start + 1;
    }

    @Override
    public int compareTo(Interval o) {
        int startComparison = Integer.compare(this.start, o.start);
        if (startComparison != 0) {
            return startComparison;
        }
        return Integer.compare(this.end, o.end); // Compare by end if starts are equal
    }
}
