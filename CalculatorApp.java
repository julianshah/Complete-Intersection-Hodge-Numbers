import javax.swing.*;
import java.awt.*;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class CalculatorApp {

    public static void main(String[] args) {
        SwingUtilities.invokeLater(CalculatorApp::createAndShowGui);
    }

    private static void createAndShowGui() {
        JFrame frame = new JFrame("My Calculator");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JTextField inputA = new JTextField(10);
        JTextField inputB = new JTextField(10);
        
        JTextField inputDimension = new JTextField(10);
        JTextField inputCodimension = new JTextField(10);
        JTextField inputDegrees = new JTextField(10);
        JButton calcButton = new JButton("Calculate");
        JLabel resultLabel = new JLabel("Result: ");

        JPanel panel = new JPanel(new GridLayout(0, 2, 8, 8));
        panel.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));
        panel.add(new JLabel("Dimension of X:"));
        panel.add(inputDimension);
        panel.add(new JLabel("Codimension of X:"));
        panel.add(inputCodimension);
        panel.add(new JLabel("Degrees:"));
        panel.add(inputDegrees);
        panel.add(calcButton);
        panel.add(resultLabel);
        
        calcButton.addActionListener(e -> {
        	try {
        		int dimension = Integer.parseInt(inputDimension.getText().trim());
        		int codimension = Integer.parseInt(inputCodimension.getText().trim());
        		List<Integer> degrees = Arrays.stream(inputDegrees.getText().split(","))
        				.map(String::trim)
        				.map(Integer::valueOf)
        				.collect(Collectors.toList());
        		if(degrees.size() != codimension)
        			throw new NumberFormatException("Not enough degrees");
        		resultLabel.setText("" + calculate(dimension, codimension, degrees));
        		frame.pack();
        	} catch (NumberFormatException ex) {
        		JOptionPane.showMessageDialog(frame, "Please enter valid numbers.",
                        "Input error", JOptionPane.ERROR_MESSAGE);
        	}
        });

        frame.getRootPane().setDefaultButton(calcButton);
        frame.add(panel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }

    /**Uses the Hilbert polynomials of Omega^p_X for each p;
     * this is enough to compute h^p,n-p since the cohomology of
     * Omega^p_X is concentrated in degrees p and n-p, and the
     * degree p part is 1-dimensional*/
    private static List<BigInteger> calculate(int dimension, int codimension, List<Integer> degrees) {
    	int ambientDimension = dimension + codimension;
    	
    	List<IntPoly> hilbertPolys = differentialHilbertPolys(dimension, codimension, degrees);
    	List<BigInteger> output = new ArrayList<BigInteger>();
    	
    	for(int p = 0; p <= dimension; p++) {
    		BigInteger euler = hilbertPolys.get(p).evaluate(0);
    		output.add(euler.subtract(BigInteger.valueOf(p%2==0 ? 1 : -1))
    				.multiply(BigInteger.valueOf((dimension-p)%2==0 ? 1 : -1)));
    	}
    	
        return output;
    }
    
    /**Computes the Hilbert polynomials of the Omega^p_X inductively
     * using exterior powers of the conormal sequence of X sitting
     * in a complete intersection one dimension higher*/
    private static List<IntPoly> differentialHilbertPolys(int dimension, int codimension, List<Integer> degrees) {
    	if(codimension == 0) {
    		List<IntPoly> output = new ArrayList<IntPoly>();
    		for(int p = 0; p <= dimension; p++) {
    			output.add(IntPoly.binomial(dimension-p).shift(dimension-p)
    					.multiply(IntPoly.binomial(p).shift(-1)));
    		}
    		return output;
    	}
    	System.out.println(degrees + "\n");
    	
    	List<Integer> choppedDegrees = new ArrayList<>(degrees);
    	Integer d = choppedDegrees.remove(codimension-1);
    	List<IntPoly> ambientHilbertPolys = differentialHilbertPolys(dimension + 1, codimension - 1, choppedDegrees);
    	
    	List<IntPoly> pullbackHilbertPolys = new ArrayList<IntPoly>();
    	for(IntPoly P : ambientHilbertPolys)
    		pullbackHilbertPolys.add(P.subtract(P.shift(-d)));
    	
    	List<IntPoly> output = new ArrayList<IntPoly>();
    	output.add(pullbackHilbertPolys.get(0));
    	for(int p = 1; p <= dimension; p++)
    		output.add(pullbackHilbertPolys.get(p).subtract(output.get(p-1).shift(-d)));
    	
    	return output;
    	
    }
}