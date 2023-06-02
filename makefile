#Makefile of different commands


###################
# Introductions
introductions:
	bash mescev/stage3/simulation_adjust_data.sh $(voc)

###################
# Projections
projections:
	bash mescev/stage3/simulation_script.sh $(voc)

###################
# run all 
all:
	bash mescev/stage3/simulation_adjust_data.sh $(voc)
	bash mescev/stage3/simulation_script.sh $(voc)

###################
# remove results 
clean:
	find results/ -type f -iname \*.csv -delete
