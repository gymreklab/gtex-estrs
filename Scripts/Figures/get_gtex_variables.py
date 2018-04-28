
def export_gtex_var():
    # Colors (use colors from GTeX e.g. http://science.sciencemag.org/content/348/6235/648.full)
    COLORS = {
        "Artery-Aorta":"salmon",
        "Artery-Tibial": "red",
        "Adipose-Subcutaneous": "darkorange",    
        "Adipose-Visceral":"orange",
        "Brain-Caudate":"lemonchiffon"   , 
        "Brain-Cerebellum":"yellow",
        "Cells-Transformedfibroblasts": "skyblue",
        "Esophagus-Mucosa": "sienna",
        "Esophagus-Muscularis":"burlywood",
        "Heart-LeftVentricle":"darkviolet",
        "Lung": "greenyellow",
        "Muscle-Skeletal": "mediumslateblue",
        "Nerve-Tibial":"gold",
        "Skin-NotSunExposed":"blue",
        "Skin-SunExposed":"cornflowerblue",
        "Thyroid":"green",
        "WholeBlood": "m",
        "permuted": "gray"
    }

    #    "Thyroid": "green",
    SHORTEN = {
        "Artery-Aorta":"Artery.A"     ,
        "Artery-Tibial": "Artery.T",
        "Adipose-Subcutaneous": "Adipose.S",    
        "Adipose-Visceral":"Adipose.V",
        "Brain-Caudate":"Caudate"   , 
        "Brain-Cerebellum":"Cerebellum",
        "Cells-Transformedfibroblasts": "Fibroblast",
        "Esophagus-Mucosa": "Mucosa",
        "Esophagus-Muscularis":"Muscularis",
        "Heart-LeftVentricle":"Heart",
        "Lung": "Lung",
        "Muscle-Skeletal": "Muscle",
        "Nerve-Tibial":"Nerve",
        "Skin-NotSunExposed": "SkinUnexposed",
        "Skin-SunExposed":"SkinLeg",
        "Thyroid":"Thyroid",
        "WholeBlood": "Blood",
        "permuted":"Permuted",
        "LCL": "LCL"
    }
    TISSUES = [item for item in list(COLORS.keys()) if item != "permuted"]
    return (COLORS, SHORTEN, TISSUES)
    
