// Macro to run ThunderSTORM on a directory containing TIF stacks/images. 
// Made by Joshua McCausland in the CJW lab, 2024
directory = getDirectory("Choose input directory");
fileList = getFileList(directory);
setBatchMode(true); // This is to run the macro in the background--do not consume extra resources to display graphics or menus.

//Define the image identifier
ident = "20230306" //can be "GFP" or whatever is unique to the actual image
idx = 0;
for(m = 0; m<fileList.length; m++) {
		if (fileList.length>0){
				file = fileList[m];
				test = indexOf(file, ident);
				if(test >= 0) {
					idx = idx + 1;
					open(file);
					run("Select All");
					
					//Camera setup parameters. Change these to match the system you are using for your experiments.
					//In our case, we only need to set pixel size. All other parameters are for estimating photons, but we just want localizations.
					run("Camera setup", "isemgain=true pixelsize=65 gainem=300.0 offset=54.0 photons2adu=4.93");
					
					//Spot detection parameters. You can optimize this by changing menu items on a test image. Once you have the settings dialed in, edit this line.
					run("Run analysis", "filter=[Wavelet filter (B-Spline)]" +
						"scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood" +
						"threshold=2*std(Wave.F1) estimator=[PSF: Gaussian] sigma=1.0 method=[Least squares]" +
						"full_image_fitting=false fitradius=3 mfaenabled=false renderer=[No Renderer]");
						
					// This adds the number iterator at the end of every image. Saves every parameter associated with the spot localization.
					if(idx == 1){
					run("Export results", "filepath=[" +directory+ "TS_results-0"+idx+".csv]"+
						"fileformat=[CSV (comma separated)] id=true frame=true sigma=true chi2=true"+
						"bkgstd=true intensity=true saveprotocol=true offset=true uncertainty=true y=true x=true");						
					}else if(idx < 10){
					run("Export results", "filepath=[" +directory+ "TS_results-0"+idx+".csv]"+
						"fileformat=[CSV (comma separated)] id=true frame=true sigma=true chi2=true"+
						"bkgstd=true intensity=true saveprotocol=false offset=true uncertainty=true y=true x=true");	
					}else{
					run("Export results", "filepath=[" +directory+ "TS_results-"+idx+".csv]"+
						"fileformat=[CSV (comma separated)] id=true frame=true sigma=true chi2=true"+
						"bkgstd=true intensity=true saveprotocol=false offset=true uncertainty=true y=true x=true");	
					}
					}
					run("Close All");
					close("Log");
				}
			}
setBatchMode("exit and display");