function slice = slice_finder(mousename, site)
mdb = initMouseDB('update', 'notext');
[~, idx] = mdb.search(mousename);
slice =  mdb.mice{idx}.phys.cell(site).file(1).Slice;