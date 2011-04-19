

plan = {
        '{maskname}':  {
                '{band}': {
                        'flats': [],
                        'arcs': [],
                        'science': [ ("{pos}-{pos}", fnA, fnB) ]
                }
        }
}

for maskname, maskplan in plan:
        for bandname, bandplan in maskplan:
                handle_flats bandplan["flats"]
                
                for position in maskplane["sience"]:
                        posname, fnA, fnB, ... = position
                        for fn in [...]:
                                Img = load_file fn 
                                fit_lambda Img
                                sky_subtract Img
                                rectify Img

                        # How to handle diffs ?
                        diff = A-B
                        residual_sky_subtract diff
                        rectify diff
