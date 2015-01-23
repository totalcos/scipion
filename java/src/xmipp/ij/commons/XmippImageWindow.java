package xmipp.ij.commons;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import java.awt.Label;
import java.awt.MenuItem;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import javax.swing.SwingUtilities;
import org.junit.experimental.theories.ParameterSignature;
import xmipp.ij.commons.XmippMenuBar.IJRequirement;
import xmipp.utils.Params;

public class XmippImageWindow extends ImageWindow implements XmippIJWindow
{

	protected XmippMenuBar menu;
	protected ImagePlusLoader ipl;
	protected Label pixelslb;
        private MaskJFrame maskfr;
        protected Params params;
    

	public XmippImageWindow(ImagePlusLoader ipl, Params params)
	{
		this(ipl, ipl.getName(), params);
	}
        
        public XmippImageWindow(ImagePlus imp, Params params)
        {
            this(imp, new XmippImageCanvas(imp), params);
        }
        
        public XmippImageWindow(ImagePlus imp, ImageCanvas canvas, Params params)
	{
		super(imp, canvas);
                this.ipl = new ImagePlusLoader(imp);
                this.params = params;
                XmippApplication.addInstance(true);
                initComponents();
	}


	public XmippImageWindow(ImagePlusLoader ipl, String title, Params params)
	{
		super(ipl.getImagePlus(), new XmippImageCanvas(ipl.getImagePlus()));
		this.ipl = ipl;
		imp.setTitle(title);
		this.params = params;
		XmippApplication.addInstance(true);
		initComponents();
                
	}
        
        protected void initComponents()
        {
            menu = new XmippMenuBar(this);
            setMenuBar(menu);
            if(XmippApplication.isScipion())
            {
                MenuItem createMaskmi = new MenuItem("Mask Manager");
                
                createMaskmi.addActionListener(new ActionListener() {

                    @Override
                    public void actionPerformed(ActionEvent ae) {
                        loadMaskFrame();
                    }
                });
                menu.advancedmn.add(createMaskmi);
                if(params.isMask())
                    loadMaskFrame();
            }
            pixelslb = new Label("                                          ");
            add(pixelslb);
            
        }

	public void openMaskToolbar()
	{
		menu.runCommand("Masks Tool Bar", new IJRequirement[] { IJRequirement.IMAGEJ });
	}

	@Override
	public void loadData()
	{
		getCanvas().loadData(this);
	}

	@Override
	public void saveDataAs(String file) throws Exception
	{
		XmippImageConverter.writeImagePlus(imp, file);
	}

	@Override
	public void saveData() throws Exception
	{
		saveDataAs(imp.getTitle());
	}

	public ImagePlusLoader getImagePlusLoader()
	{
		return ipl;
	}

	@Override
	public boolean isVolume()
	{
		return false;
	}

	@Override
	public boolean isStack()
	{
		return false;
	}

	//overwriting ImageJ event to avoid switching menu
	public void windowActivated(WindowEvent e)
	{
		//		if (IJ.isMacintosh())
		//			this.setMenuBar(Menus.getMenuBar());
		if (IJ.debugMode)
			IJ.write(imp.getTitle() + ": Activated");
		if (!closed)
		{
			//ic.requestFocus();
			WindowManager.setCurrentWindow(this);
		}
	}

	public XmippImageCanvas getCanvas()
	{
		return ((XmippImageCanvas) super.getCanvas());
	}

	public void showPixels(int x, int y, int[] pixels)
	{
                String text;
		String value = "";
		switch (imp.getType())
		{
		case ImagePlus.GRAY8:
		case ImagePlus.GRAY16:
			double cValue = imp.getCalibration().getCValue(pixels[0]);
			if (cValue == pixels[0])
				value = String.valueOf(pixels[0]);
			else
				value = IJ.d2s(cValue) + " (" + pixels[0] + ")";
                        text = String.format("x=%s, y=%s, value=%-5s", x, y, value);
                        pixelslb.setText(text);
			break;
		case ImagePlus.GRAY32:
                        text = String.format("x=%s, y=%s, value=%.2f", x, y, Float.intBitsToFloat(pixels[0]));
                        pixelslb.setText(text);
			break;
		case ImagePlus.COLOR_256:
		case ImagePlus.COLOR_RGB:
                        value =  pixels[0] + "," + pixels[1] + "," + pixels[2];
                        text = String.format("x=%s, y=%s, value=%-15s", x, y, value);
                        pixelslb.setText(text);
			break;
		
		}
		
	}
        
        
        
        @Override
	public void windowClosing(WindowEvent e) {
            
            super.windowClosing(e);
            XmippApplication.removeInstance(true);
            if(maskfr != null)
            {
                maskfr.setVisible(false);
                maskfr.dispose();
                
            }
		
	}
        
       
        public Params getParams()
        {
            return params;
        }

        public void setParams(Params params)
        {
            this.params = params;
        }

    
        protected void loadMaskFrame()
        {
            if(maskfr == null || !maskfr.isVisible())
            {
                SwingUtilities.invokeLater(new Runnable() {

                    @Override
                    public void run() {
                        maskfr = new MaskJFrame(XmippImageWindow.this);
                    }
                });
                
                imp.getCanvas().addMouseListener(new MouseListener() {

                    @Override
                    public void mouseClicked(MouseEvent me) {
                    }

                    @Override
                    public void mousePressed(MouseEvent me) {
                    }

                    @Override
                    public void mouseReleased(MouseEvent me) {
                        maskfr.refreshMask();
                    }

                    @Override
                    public void mouseEntered(MouseEvent me) {
                    }

                    @Override
                    public void mouseExited(MouseEvent me) {
                    }
                });
                imp.getCanvas().addComponentListener(new java.awt.event.ComponentAdapter()
		{
			public void componentResized(ComponentEvent e)
			{
				maskfr.refreshMask();
			}
		});
            }
        }

}// class XmippImageWindow
