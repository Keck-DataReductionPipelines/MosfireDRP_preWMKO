

plan = {
    '{maskname}':  {
        '{band}': {
            'flats': [],
            'arcs': [],
            'science': [ ("{pos}-{pos}", fnA, fnB) ]
        }
    }
}


high_level function ------
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


def fit_lambda img, hdr:
    
    solutions = []
    for each slit:
        y0 = central_pixel(slit)
        [order, y0, alpha, sinbeta, gamma, delta] = 
            guess_wavelength_solution(img, slit)
        sol_1d = fit_wavelength(data, y0, [order, y0, alpha, 
            sinbeta, gamma, delta], slit)
        sol_2d = fit_outwards(data, y0, sol1d, slit)

        solutions.append 

    return solutions

