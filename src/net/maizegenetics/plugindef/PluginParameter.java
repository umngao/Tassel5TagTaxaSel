package net.maizegenetics.plugindef;

import com.google.common.collect.Range;

import static net.maizegenetics.plugindef.AbstractPlugin.convert;

import org.apache.log4j.Logger;

/**
 * Defines the attributes of parameters to be used in the plugins
 *
 * @param <T>
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 *
 */
public final class PluginParameter<T> {

    private static final Logger myLogger = Logger.getLogger(PluginParameter.class);

    private final String myGuiName;
    private final String myUnits;
    private final String myCmdLineName;
    private final String myDescription;
    private final Range<Comparable<T>> myRange;
    private final T myDefaultValue;
    private final T myValue;
    private final boolean myRequired;
    private final Class<T> myClass;
    private final PluginParameter<?> myDependentOnParameter;
    private final Object myDependentOnParameterValue;

    public enum PARAMETER_TYPE {

        NA, IN_FILE, OUT_FILE, IN_DIR, OUT_DIR, GENOTYPE_TABLE, TAXA_LIST, POSITION_LIST
    };
    private final PARAMETER_TYPE myParameterType;

    private PluginParameter(String guiName, String guiUnits, String cmdLineName,
            String description, Range<Comparable<T>> range, T defaultValue, T value, boolean required, PARAMETER_TYPE fileType, PluginParameter<?> dependentOnParameter, Object dependentOnParameterValue, Class<T> type) {
        myGuiName = guiName;
        myUnits = guiUnits;
        myCmdLineName = cmdLineName;
        myDescription = description;
        myRange = range;
        myDefaultValue = defaultValue;
        if (value == null) {
            myValue = defaultValue;
        } else {
            myValue = value;
        }

        if ((myRange != null) && (myValue != null)) {
            if (!acceptsValue((Comparable<T>) myValue)) {
                StringBuilder builder = new StringBuilder();
                builder.append("PluginParameter: init: " + myCmdLineName + " value: " + value.toString() + " outside range: ");
                if (valueType().isEnum()) {
                    builder.append(" [");
                    T[] values = valueType().getEnumConstants();
                    for (int i = 0; i < values.length; i++) {
                        if (i != 0) {
                            builder.append(" ");
                        }
                        builder.append(values[i].toString());
                    }
                    builder.append("]");
                } else {
                    builder.append(myRange.toString());
                }
                throw new IllegalArgumentException(builder.toString());
            }
        }

        myRequired = required;
        if ((myDefaultValue != null) && (myRequired)) {
            throw new IllegalArgumentException("PluginParameter: init: " + myCmdLineName + " shouldn't have default value and be required.");
        }
        myClass = type;
        myParameterType = fileType;
        myDependentOnParameter = dependentOnParameter;
        myDependentOnParameterValue = dependentOnParameterValue;
    }

    /**
     * Use these to change the value of an existing parameter, e.g. after a user
     * changes the value. Otherwise use the Builder to create the parameter
     *
     * @param oldParameter
     * @param newValue
     */
    public PluginParameter(PluginParameter<T> oldParameter, T newValue) {
        this(oldParameter.myGuiName, oldParameter.myUnits, oldParameter.myCmdLineName,
                oldParameter.myDescription, oldParameter.myRange, oldParameter.myDefaultValue, newValue,
                oldParameter.myRequired, oldParameter.myParameterType, oldParameter.dependentOnParameter(),
                oldParameter.dependentOnParameterValue(), oldParameter.myClass);
    }

    public String guiName() {
        return myGuiName;
    }

    public String units() {
        return myUnits;
    }

    public String cmdLineName() {
        return myCmdLineName;
    }

    public String description() {
        return myDescription;
    }

    public Range<Comparable<T>> range() {
        return myRange;
    }

    public boolean acceptsValue(Object value) {
        try {
            return (myRange == null) || (myRange.contains((Comparable<T>) value));
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            return false;
        }
    }

    public boolean acceptsValue(String input) {
        Comparable<T> value = (Comparable<T>) convert(input, valueType());
        return (myRange == null) || (myRange.contains(value));
    }

    public T value() {
        return myValue;
    }

    public T defaultValue() {
        return myDefaultValue;
    }

    public boolean required() {
        return myRequired;
    }

    public Class<T> valueType() {
        return myClass;
    }

    public PARAMETER_TYPE parameterType() {
        return myParameterType;
    }

    public PluginParameter<?> dependentOnParameter() {
        return myDependentOnParameter;
    }

    public Object dependentOnParameterValue() {
        return myDependentOnParameterValue;
    }

    public boolean isEmpty() {
        return (myValue == null) || (myValue.toString().trim().length() == 0);
    }

    public static class Builder<T> {

        private String myGuiName;
        private String myUnits = "";
        private final String myCmdLineName;
        private String myDescription = "";
        private Range<Comparable<T>> myRange = null;
        private final T myDefaultValue;
        private boolean myIsRequired = false;
        private final Class<T> myClass;
        private PARAMETER_TYPE myParameterType = PARAMETER_TYPE.NA;
        private PluginParameter<?> myDependentOnParameter = null;
        private Object myDependentOnParameterValue = null;

        public Builder(String cmdLineName, T defaultValue, Class<T> type) {
            myCmdLineName = cmdLineName;
            myDefaultValue = defaultValue;
            myClass = type;
        }

        public Builder<T> units(String units) {
            myUnits = units;
            return this;
        }

        public Builder<T> description(String description) {
            myDescription = description;
            return this;
        }

        public Builder<T> range(Range<Comparable<T>> range) {
            myRange = range;
            return this;
        }

        public Builder<T> required(boolean required) {
            myIsRequired = required;
            return this;
        }

        public Builder<T> guiName(String guiName) {
            myGuiName = guiName;
            return this;
        }

        public Builder<T> inFile() {
            myParameterType = PARAMETER_TYPE.IN_FILE;
            return this;
        }

        public Builder<T> outFile() {
            myParameterType = PARAMETER_TYPE.OUT_FILE;
            return this;
        }

        public Builder<T> inDir() {
            myParameterType = PARAMETER_TYPE.IN_DIR;
            return this;
        }

        public Builder<T> outDir() {
            myParameterType = PARAMETER_TYPE.OUT_DIR;
            return this;
        }

        public Builder<T> genotypeTable() {
            myParameterType = PARAMETER_TYPE.GENOTYPE_TABLE;
            return this;
        }

        public Builder<T> taxaList() {
            myParameterType = PARAMETER_TYPE.TAXA_LIST;
            return this;
        }

        public Builder<T> positionList() {
            myParameterType = PARAMETER_TYPE.POSITION_LIST;
            return this;
        }

        public Builder<T> dependentOnParameter(PluginParameter<?> parameter) {
            if (Boolean.class.isAssignableFrom(parameter.valueType())) {
                return dependentOnParameter(parameter, true);
            } else {
                throw new IllegalArgumentException("PluginParameter: dependentOnParameter: no default value for: " + parameter.valueType().getName());
            }
        }

        public Builder<T> dependentOnParameter(PluginParameter<?> parameter, Object value) {
            myDependentOnParameter = parameter;
            myDependentOnParameterValue = value;
            return this;
        }

        public PluginParameter<T> build() {
            if ((myGuiName == null) || (myGuiName.isEmpty())) {
                StringBuilder builder = new StringBuilder();
                builder.append(Character.toUpperCase(myCmdLineName.charAt(0)));
                for (int i = 1; i < myCmdLineName.length(); i++) {
                    char current = myCmdLineName.charAt(i);
                    if (Character.isUpperCase(current)) {
                        builder.append(" ");
                    }
                    builder.append(current);
                }
                myGuiName = builder.toString();
            }
            if (myDescription.isEmpty()) {
                myDescription = myGuiName;
            }
            return new PluginParameter<>(myGuiName, myUnits, myCmdLineName,
                    myDescription, myRange, myDefaultValue, null, myIsRequired,
                    myParameterType, myDependentOnParameter,
                    myDependentOnParameterValue, myClass);
        }
    }
}
