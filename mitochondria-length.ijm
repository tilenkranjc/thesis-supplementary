
dir1 = getDirectory("Choose Source Directory ");
format = getFormat();
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);
setBatchMode(true);
for (i=0; i<list.length; i++) {
    showProgress(i+1, list.length);
    open(dir1+list[i]);
    if (format=="8-bit TIFF" || format=="GIF")
        convertTo8Bit();
    // Rolling ball with 1 pixel
    run("Subtract Background...", "rolling=1");
    // Automatic threshold
    setAutoThreshold("Otsu dark");
    setOption("BlackBackground", false);
    run("Convert to Mask");
    // Reset the scale of the image to pixels
    run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
    setOption("BlackBackground", false);
    // Skeletonize the binary mask
    run("Skeletonize");
    // Measure length of the skeletons that are longer than 7 px
    run("Analyze Particles...", "size=7-300 display clear");
    // Save results file
    saveAs("Results", dir2+list[i]+".tab");
    close();
}

function getFormat() {
    formats = newArray("TIFF", "8-bit TIFF", "JPEG", "GIF", "PNG",
        "PGM", "BMP", "FITS", "Text Image", "ZIP", "Raw");
    Dialog.create("Batch Convert");
    Dialog.addChoice("Convert to: ", formats, "TIFF");
    Dialog.show();
    return Dialog.getChoice();
}

function convertTo8Bit() {
    if (bitDepth==24)
        run("8-bit Color", "number=256");
    else
        run("8-bit");
}