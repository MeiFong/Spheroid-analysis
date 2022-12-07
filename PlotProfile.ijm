// ask user to select a folder
dir = getDirectory("Select A folder");
// get the list of files (& folders) in it
fileList = getFileList(dir);
// prepare a folder to output the images
output_dir = dir + File.separator + "output" + File.separator ;
File.makeDirectory(output_dir);


//activate batch mode
setBatchMode(true);

// LOOP to process the list of files
for (i = 0; i < lengthOf(fileList); i++) {
	// define the "path" 
	// by concatenation of dir and the i element of the array fileList
	current_imagePath = dir+fileList[i];
	// check that the currentFile is not a directory
	if (!File.isDirectory(current_imagePath)){

		open(current_imagePath);
		makeRectangle(0, 0, 1000, 1000);
		run("Clear Results");
  		
  		run("Plot Profile");
		Plot.showValues();
  		
  		// Save as spreadsheet compatible text file
  		current_imagename = getTitle();
  		name = current_imagename + ".csv";
  		saveAs("Results", output_dir + name);

		
	}	
		// make sure to close every images befores opening the next one
		run("Close All");
}
setBatchMode(false);