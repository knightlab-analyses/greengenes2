#!/bin/bash
#SBATCH -J effect
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 64g
#SBATCH --time 4:00:00

source activate qiime2-2022.2

set -e


qiime evident beta-effect-size-by-category \
    --i-beta-diversity ${dm1} \
    --m-sample-metadata-file ${metadata} \
    --p-columns age_cat \
    --p-columns thdmi_cohort \
    --p-columns types_of_plants \
    --p-columns acid_reflux \
    --p-columns alcohol_consumption \
    --p-columns alcohol_frequency \
    --p-columns alcohol_types_beercider \
    --p-columns alcohol_types_red_wine \
    --p-columns alcohol_types_spiritshard_alcohol \
    --p-columns alcohol_types_white_wine \
    --p-columns allergic_to_i_have_no_food_allergies_that_i_know_of \
    --p-columns artificial_sweeteners \
    --p-columns beet_frequency \
    --p-columns bmi_cat \
    --p-columns bowel_movement_frequency \
    --p-columns cat \
    --p-columns clinical_condition \
    --p-columns consume_animal_products_abx \
    --p-columns cosmetics_frequency \
    --p-columns covid_anxious \
    --p-columns covid_chronic_conditions_seasonal_allergies \
    --p-columns covid_depressed \
    --p-columns covid_difficulty_falling_asleep \
    --p-columns covid_difficulty_staying_asleep \
    --p-columns covid_family_members_none_of_the_above \
    --p-columns covid_family_members_put_into_self_quarantine_without_symptoms_eg_due_to_possible_exposure \
    --p-columns covid_happened_to_you_none_of_the_above \
    --p-columns covid_interest_pleasure \
    --p-columns covid_left_home \
    --p-columns covid_level_of_wellbeing \
    --p-columns covid_likely_exposure_no \
    --p-columns covid_quality_of_life \
    --p-columns covid_sleep_interference \
    --p-columns covid_sleep_pattern \
    --p-columns covid_symptoms_fatigue \
    --p-columns covid_waking_up_early \
    --p-columns covid_worried_sleep \
    --p-columns covid_worrying \
    --p-columns deodorant_use \
    --p-columns dog \
    --p-columns drinking_water_source \
    --p-columns drinks_per_session \
    --p-columns exercise_frequency \
    --p-columns exercise_location \
    --p-columns fed_as_infant \
    --p-columns fermented_plant_frequency \
    --p-columns flossing_frequency \
    --p-columns flu_vaccine_date \
    --p-columns fruit_frequency \
    --p-columns high_fat_red_meat_frequency \
    --p-columns homecooked_meals_frequency \
    --p-columns ibs \
    --p-columns level_of_education \
    --p-columns livingwith \
    --p-columns meat_eggs_frequency \
    --p-columns milk_cheese_frequency \
    --p-columns milk_substitute_frequency \
    --p-columns multivitamin \
    --p-columns nail_biter \
    --p-columns non_food_allergies_unspecified \
    --p-columns olive_oil \
    --p-columns one_liter_of_water_a_day_frequency \
    --p-columns other_supplement_frequency \
    --p-columns plant_protein_frequency \
    --p-columns pool_frequency \
    --p-columns poultry_frequency \
    --p-columns prepared_meals_frequency \
    --p-columns probiotic_frequency \
    --p-columns red_meat_frequency \
    --p-columns salted_snacks_frequency \
    --p-columns seafood_frequency \
    --p-columns seasonal_allergies \
    --p-columns sleep_duration \
    --p-columns softener \
    --p-columns specialized_diet_i_do_not_eat_a_specialized_diet \
    --p-columns sugar_sweetened_drink_frequency \
    --p-columns sugary_sweets_frequency \
    --p-columns tonsils_removed \
    --p-columns vegetable_frequency \
    --p-columns vitamin_b_supplement_frequency \
    --p-columns vitamin_d_supplement_frequency \
    --p-columns vivid_dreams \
    --p-columns whole_eggs \
    --p-columns whole_grain_frequency \
    --o-effect-size-results ${output}.dm1.qza \
    --verbose

qiime evident beta-effect-size-by-category \
    --i-beta-diversity ${dm2} \
    --m-sample-metadata-file ${metadata} \
    --p-columns age_cat \
    --p-columns thdmi_cohort \
    --p-columns types_of_plants \
    --p-columns acid_reflux \
    --p-columns alcohol_consumption \
    --p-columns alcohol_frequency \
    --p-columns alcohol_types_beercider \
    --p-columns alcohol_types_red_wine \
    --p-columns alcohol_types_spiritshard_alcohol \
    --p-columns alcohol_types_white_wine \
    --p-columns allergic_to_i_have_no_food_allergies_that_i_know_of \
    --p-columns artificial_sweeteners \
    --p-columns beet_frequency \
    --p-columns bmi_cat \
    --p-columns bowel_movement_frequency \
    --p-columns cat \
    --p-columns clinical_condition \
    --p-columns consume_animal_products_abx \
    --p-columns cosmetics_frequency \
    --p-columns covid_anxious \
    --p-columns covid_chronic_conditions_seasonal_allergies \
    --p-columns covid_depressed \
    --p-columns covid_difficulty_falling_asleep \
    --p-columns covid_difficulty_staying_asleep \
    --p-columns covid_family_members_none_of_the_above \
    --p-columns covid_family_members_put_into_self_quarantine_without_symptoms_eg_due_to_possible_exposure \
    --p-columns covid_happened_to_you_none_of_the_above \
    --p-columns covid_interest_pleasure \
    --p-columns covid_left_home \
    --p-columns covid_level_of_wellbeing \
    --p-columns covid_likely_exposure_no \
    --p-columns covid_quality_of_life \
    --p-columns covid_sleep_interference \
    --p-columns covid_sleep_pattern \
    --p-columns covid_symptoms_fatigue \
    --p-columns covid_waking_up_early \
    --p-columns covid_worried_sleep \
    --p-columns covid_worrying \
    --p-columns deodorant_use \
    --p-columns dog \
    --p-columns drinking_water_source \
    --p-columns drinks_per_session \
    --p-columns exercise_frequency \
    --p-columns exercise_location \
    --p-columns fed_as_infant \
    --p-columns fermented_plant_frequency \
    --p-columns flossing_frequency \
    --p-columns flu_vaccine_date \
    --p-columns fruit_frequency \
    --p-columns high_fat_red_meat_frequency \
    --p-columns homecooked_meals_frequency \
    --p-columns ibs \
    --p-columns level_of_education \
    --p-columns livingwith \
    --p-columns meat_eggs_frequency \
    --p-columns milk_cheese_frequency \
    --p-columns milk_substitute_frequency \
    --p-columns multivitamin \
    --p-columns nail_biter \
    --p-columns non_food_allergies_unspecified \
    --p-columns olive_oil \
    --p-columns one_liter_of_water_a_day_frequency \
    --p-columns other_supplement_frequency \
    --p-columns plant_protein_frequency \
    --p-columns pool_frequency \
    --p-columns poultry_frequency \
    --p-columns prepared_meals_frequency \
    --p-columns probiotic_frequency \
    --p-columns red_meat_frequency \
    --p-columns salted_snacks_frequency \
    --p-columns seafood_frequency \
    --p-columns seasonal_allergies \
    --p-columns sleep_duration \
    --p-columns softener \
    --p-columns specialized_diet_i_do_not_eat_a_specialized_diet \
    --p-columns sugar_sweetened_drink_frequency \
    --p-columns sugary_sweets_frequency \
    --p-columns tonsils_removed \
    --p-columns vegetable_frequency \
    --p-columns vitamin_b_supplement_frequency \
    --p-columns vitamin_d_supplement_frequency \
    --p-columns vivid_dreams \
    --p-columns whole_eggs \
    --p-columns whole_grain_frequency \
    --o-effect-size-results ${output}.dm2.qza \
    --verbose

qiime evident visualize-results \
    --i-results ${output}.dm1.qza \
    --o-visualization ${output}.dm1.qzv

qiime evident visualize-results \
    --i-results ${output}.dm2.qza \
    --o-visualization ${output}.dm2.qzv

python correlate_effect_sizes.py \
    --dm1 ${output}.dm1.qza \
    --dm2 ${output}.dm2.qza \
    --dm1label ${dm1label} \
    --dm2label ${dm2label} \
    --output ${output}.scatter.pdf
