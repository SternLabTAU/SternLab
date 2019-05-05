#! /usr/local/python_anaconda/bin/python3.4



BALTIMORE = {"Arena":"ssRNA(-)",
             "Calici":"ssRNA(+)",
             "Corona":"ssRNA(+)",
             "Filo":"ssRNA(-)",
             "Flavi":"ssRNA(+)",
             "Hanta":"ssRNA(-)",
             "Hepe":"ssRNA(+)",
             "Herpes":"dsDNA",
             "HIV1":"RT",
             "Nairo":"ssRNA(-)",
             "Paramyxo":"ssRNA(-)",
             "Peribunya":"ssRNA(-)",
             "Phenui":"ssRNA(-)",
             "Orthomyxo":"ssRNA(-)",
             "Picorna":"ssRNA(+)",
             "Pox":"dsDNA",
             "Reo":"dsRNA",
             "Rhabdo":"ssRNA(-)",
             "Toga":"ssRNA(+)"}

def get_baltimore_classifiaction(family):
    return BALTIMORE[family]