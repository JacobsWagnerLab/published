%% Select a color and get its normalized RGB values.
%
% -About-
%   One often needs to select a color for various plotting. However, the
%   built-in MATLAB color palettes are limited. Based on hundreds of named
%   colors, we can extract their RGB values for better looking plots.
%
% -Input-
%   - color: a char array of a color name, note the space in the name must
%            be replaced with '_' (underscores).
%
% -Output-
%   - rgb:   normalized RGB values of the input color. If the color name
%            was not found, RGB of "black" is output
%
% -Example-
%   - myParticle.selectColor('yale_blue')
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function rgb = selectColor(~,color)
colors = ...
   {'air_force_blue_raf'
    'air_force_blue_usaf'
    'air_superiority_blue'
    'alabama_crimson'
    'alice_blue'
    'alizarin_crimson'
    'alloy_orange'
    'almond'
    'amaranth'
    'amber'
    'amber_sae_ece'
    'american_rose'
    'amethyst'
    'android_green'
    'anti_flash_white'
    'antique_brass'
    'antique_fuchsia'
    'antique_ruby'
    'antique_white'
    'ao_english'
    'apple_green'
    'apricot'
    'aqua'
    'aquamarine'
    'army_green'
    'arsenic'
    'arylide_yellow'
    'ash_grey'
    'asparagus'
    'atomic_tangerine'
    'auburn'
    'aureolin'
    'aurometalsaurus'
    'avocado'
    'azure'
    'azure_mist_web'
    'baby_blue'
    'baby_blue_eyes'
    'baby_pink'
    'ball_blue'
    'banana_mania'
    'banana_yellow'
    'barn_red'
    'battleship_grey'
    'bazaar'
    'beau_blue'
    'beaver'
    'beige'
    'big_dip_o_ruby'
    'bisque'
    'bistre'
    'bittersweet'
    'bittersweet_shimmer'
    'black'
    'black_bean'
    'black_leather_jacket'
    'black_olive'
    'blanched_almond'
    'blast_off_bronze'
    'bleu_de_france'
    'blizzard_blue'
    'blond'
    'blue'
    'blue_bell'
    'blue_crayola'
    'blue_gray'
    'blue_green'
    'blue_munsell'
    'blue_ncs'
    'blue_pigment'
    'blue_ryb'
    'blue_sapphire'
    'blue_violet'
    'blush'
    'bole'
    'bondi_blue'
    'bone'
    'boston_university_red'
    'bottle_green'
    'boysenberry'
    'brandeis_blue'
    'brass'
    'brick_red'
    'bright_cerulean'
    'bright_green'
    'bright_lavender'
    'bright_maroon'
    'bright_pink'
    'bright_turquoise'
    'bright_ube'
    'brilliant_lavender'
    'brilliant_rose'
    'brink_pink'
    'british_racing_green'
    'bronze'
    'brown_traditional'
    'brown_web'
    'bubble_gum'
    'bubbles'
    'buff'
    'bulgarian_rose'
    'burgundy'
    'burlywood'
    'burnt_orange'
    'burnt_sienna'
    'burnt_umber'
    'byzantine'
    'byzantium'
    'cadet'
    'cadet_blue'
    'cadet_grey'
    'cadmium_green'
    'cadmium_orange'
    'cadmium_red'
    'cadmium_yellow'
    'caf_au_lait'
    'caf_noir'
    'cal_poly_green'
    'cambridge_blue'
    'camel'
    'cameo_pink'
    'camouflage_green'
    'canary_yellow'
    'candy_apple_red'
    'candy_pink'
    'capri'
    'caput_mortuum'
    'cardinal'
    'caribbean_green'
    'carmine'
    'carmine_m_p'
    'carmine_pink'
    'carmine_red'
    'carnation_pink'
    'carnelian'
    'carolina_blue'
    'carrot_orange'
    'catalina_blue'
    'ceil'
    'celadon'
    'celadon_blue'
    'celadon_green'
    'celeste_colour'
    'celestial_blue'
    'cerise'
    'cerise_pink'
    'cerulean'
    'cerulean_blue'
    'cerulean_frost'
    'cg_blue'
    'cg_red'
    'chamoisee'
    'champagne'
    'charcoal'
    'charm_pink'
    'chartreuse_traditional'
    'chartreuse_web'
    'cherry'
    'cherry_blossom_pink'
    'chestnut'
    'china_pink'
    'china_rose'
    'chinese_red'
    'chocolate_traditional'
    'chocolate_web'
    'chrome_yellow'
    'cinereous'
    'cinnabar'
    'cinnamon'
    'citrine'
    'classic_rose'
    'cobalt'
    'cocoa_brown'
    'coffee'
    'columbia_blue'
    'congo_pink'
    'cool_black'
    'cool_grey'
    'copper'
    'copper_crayola'
    'copper_penny'
    'copper_red'
    'copper_rose'
    'coquelicot'
    'coral'
    'coral_pink'
    'coral_red'
    'cordovan'
    'corn'
    'cornell_red'
    'cornflower_blue'
    'cornsilk'
    'cosmic_latte'
    'cotton_candy'
    'cream'
    'crimson'
    'crimson_glory'
    'cyan'
    'cyan_process'
    'daffodil'
    'dandelion'
    'dark_blue'
    'dark_brown'
    'dark_byzantium'
    'dark_candy_apple_red'
    'dark_cerulean'
    'dark_chestnut'
    'dark_coral'
    'dark_cyan'
    'dark_electric_blue'
    'dark_goldenrod'
    'dark_gray'
    'dark_green'
    'dark_imperial_blue'
    'dark_jungle_green'
    'dark_khaki'
    'dark_lava'
    'dark_lavender'
    'dark_magenta'
    'dark_midnight_blue'
    'dark_olive_green'
    'dark_orange'
    'dark_orchid'
    'dark_pastel_blue'
    'dark_pastel_green'
    'dark_pastel_purple'
    'dark_pastel_red'
    'dark_pink'
    'dark_powder_blue'
    'dark_raspberry'
    'dark_red'
    'dark_salmon'
    'dark_scarlet'
    'dark_sea_green'
    'dark_sienna'
    'dark_slate_blue'
    'dark_slate_gray'
    'dark_spring_green'
    'dark_tan'
    'dark_tangerine'
    'dark_taupe'
    'dark_terra_cotta'
    'dark_turquoise'
    'dark_violet'
    'dark_yellow'
    'dartmouth_green'
    'davy_s_grey'
    'debian_red'
    'deep_carmine'
    'deep_carmine_pink'
    'deep_carrot_orange'
    'deep_cerise'
    'deep_champagne'
    'deep_chestnut'
    'deep_coffee'
    'deep_fuchsia'
    'deep_jungle_green'
    'deep_lilac'
    'deep_magenta'
    'deep_peach'
    'deep_pink'
    'deep_ruby'
    'deep_saffron'
    'deep_sky_blue'
    'deep_tuscan_red'
    'denim'
    'desert'
    'desert_sand'
    'dim_gray'
    'dodger_blue'
    'dogwood_rose'
    'dollar_bill'
    'drab'
    'duke_blue'
    'earth_yellow'
    'ebony'
    'ecru'
    'eggplant'
    'eggshell'
    'egyptian_blue'
    'electric_blue'
    'electric_crimson'
    'electric_cyan'
    'electric_green'
    'electric_indigo'
    'electric_lavender'
    'electric_lime'
    'electric_purple'
    'electric_ultramarine'
    'electric_violet'
    'electric_yellow'
    'emerald'
    'english_lavender'
    'eton_blue'
    'fallow'
    'falu_red'
    'fandango'
    'fashion_fuchsia'
    'fawn'
    'feldgrau'
    'fern_green'
    'ferrari_red'
    'field_drab'
    'fire_engine_red'
    'firebrick'
    'flame'
    'flamingo_pink'
    'flavescent'
    'flax'
    'floral_white'
    'fluorescent_orange'
    'fluorescent_pink'
    'fluorescent_yellow'
    'folly'
    'forest_green_traditional'
    'forest_green_web'
    'french_beige'
    'french_blue'
    'french_lilac'
    'french_lime'
    'french_raspberry'
    'french_rose'
    'fuchsia'
    'fuchsia_crayola'
    'fuchsia_pink'
    'fuchsia_rose'
    'fulvous'
    'fuzzy_wuzzy'
    'gainsboro'
    'gamboge'
    'ghost_white'
    'ginger'
    'glaucous'
    'glitter'
    'gold_metallic'
    'gold_web_golden'
    'golden_brown'
    'golden_poppy'
    'golden_yellow'
    'goldenrod'
    'granny_smith_apple'
    'gray'
    'gray_asparagus'
    'gray_html_css_gray'
    'gray_x11_gray'
    'green_color_wheel_x11_green'
    'green_crayola'
    'green_html_css_green'
    'green_munsell'
    'green_ncs'
    'green_pigment'
    'green_ryb'
    'green_yellow'
    'grullo'
    'guppie_green'
    'halay_be'
    'han_blue'
    'han_purple'
    'hansa_yellow'
    'harlequin'
    'harvard_crimson'
    'harvest_gold'
    'heart_gold'
    'heliotrope'
    'hollywood_cerise'
    'honeydew'
    'honolulu_blue'
    'hooker_s_green'
    'hot_magenta'
    'hot_pink'
    'hunter_green'
    'iceberg'
    'icterine'
    'imperial_blue'
    'inchworm'
    'india_green'
    'indian_red'
    'indian_yellow'
    'indigo'
    'indigo_dye'
    'indigo_web'
    'international_klein_blue'
    'international_orange_aerospace'
    'international_orange_engineering'
    'international_orange_golden_gate_bridge'
    'iris'
    'isabelline'
    'islamic_green'
    'ivory'
    'jade'
    'jasmine'
    'jasper'
    'jazzberry_jam'
    'jet'
    'jonquil'
    'june_bud'
    'jungle_green'
    'kelly_green'
    'kenyan_copper'
    'khaki_html_css_khaki'
    'khaki_x11_light_khaki'
    'ku_crimson'
    'la_salle_green'
    'languid_lavender'
    'lapis_lazuli'
    'laser_lemon'
    'laurel_green'
    'lava'
    'lavender_blue'
    'lavender_blush'
    'lavender_floral'
    'lavender_gray'
    'lavender_indigo'
    'lavender_magenta'
    'lavender_mist'
    'lavender_pink'
    'lavender_purple'
    'lavender_rose'
    'lavender_web'
    'lawn_green'
    'lemon'
    'lemon_chiffon'
    'lemon_lime'
    'licorice'
    'light_apricot'
    'light_blue'
    'light_brown'
    'light_carmine_pink'
    'light_coral'
    'light_cornflower_blue'
    'light_crimson'
    'light_cyan'
    'light_fuchsia_pink'
    'light_goldenrod_yellow'
    'light_gray'
    'light_green'
    'light_khaki'
    'light_pastel_purple'
    'light_pink'
    'light_red_ochre'
    'light_salmon'
    'light_salmon_pink'
    'light_sea_green'
    'light_sky_blue'
    'light_slate_gray'
    'light_taupe'
    'light_thulian_pink'
    'light_yellow'
    'lilac'
    'lime_color_wheel'
    'lime_green'
    'lime_web_x11_green'
    'limerick'
    'lincoln_green'
    'linen'
    'lion'
    'little_boy_blue'
    'liver'
    'lust'
    'magenta'
    'magenta_dye'
    'magenta_process'
    'magic_mint'
    'magnolia'
    'mahogany'
    'maize'
    'majorelle_blue'
    'malachite'
    'manatee'
    'mango_tango'
    'mantis'
    'mardi_gras'
    'maroon_crayola'
    'maroon_html_css'
    'maroon_x11'
    'mauve'
    'mauve_taupe'
    'mauvelous'
    'maya_blue'
    'meat_brown'
    'medium_aquamarine'
    'medium_blue'
    'medium_candy_apple_red'
    'medium_carmine'
    'medium_champagne'
    'medium_electric_blue'
    'medium_jungle_green'
    'medium_lavender_magenta'
    'medium_orchid'
    'medium_persian_blue'
    'medium_purple'
    'medium_red_violet'
    'medium_ruby'
    'medium_sea_green'
    'medium_slate_blue'
    'medium_spring_bud'
    'medium_spring_green'
    'medium_taupe'
    'medium_turquoise'
    'medium_tuscan_red'
    'medium_vermilion'
    'medium_violet_red'
    'mellow_apricot'
    'mellow_yellow'
    'melon'
    'midnight_blue'
    'midnight_green_eagle_green'
    'mikado_yellow'
    'mint'
    'mint_cream'
    'mint_green'
    'misty_rose'
    'moccasin'
    'mode_beige'
    'moonstone_blue'
    'mordant_red_19'
    'moss_green'
    'mountain_meadow'
    'mountbatten_pink'
    'msu_green'
    'mulberry'
    'mustard'
    'myrtle'
    'nadeshiko_pink'
    'napier_green'
    'naples_yellow'
    'navajo_white'
    'navy_blue'
    'neon_carrot'
    'neon_fuchsia'
    'neon_green'
    'new_york_pink'
    'non_photo_blue'
    'north_texas_green'
    'ocean_boat_blue'
    'ochre'
    'office_green'
    'old_gold'
    'old_lace'
    'old_lavender'
    'old_mauve'
    'old_rose'
    'olive'
    'olive_drab_7'
    'olive_drab_web_olive_drab_3'
    'olivine'
    'onyx'
    'opera_mauve'
    'orange_color_wheel'
    'orange_peel'
    'orange_red'
    'orange_ryb'
    'orange_web_color'
    'orchid'
    'otter_brown'
    'ou_crimson_red'
    'outer_space'
    'outrageous_orange'
    'oxford_blue'
    'pakistan_green'
    'palatinate_blue'
    'palatinate_purple'
    'pale_aqua'
    'pale_blue'
    'pale_brown'
    'pale_carmine'
    'pale_cerulean'
    'pale_chestnut'
    'pale_copper'
    'pale_cornflower_blue'
    'pale_gold'
    'pale_goldenrod'
    'pale_green'
    'pale_lavender'
    'pale_magenta'
    'pale_pink'
    'pale_plum'
    'pale_red_violet'
    'pale_robin_egg_blue'
    'pale_silver'
    'pale_spring_bud'
    'pale_taupe'
    'pale_violet_red'
    'pansy_purple'
    'papaya_whip'
    'paris_green'
    'pastel_blue'
    'pastel_brown'
    'pastel_gray'
    'pastel_green'
    'pastel_magenta'
    'pastel_orange'
    'pastel_pink'
    'pastel_purple'
    'pastel_red'
    'pastel_violet'
    'pastel_yellow'
    'patriarch'
    'payne_s_grey'
    'peach'
    'peach_crayola'
    'peach_orange'
    'peach_puff'
    'peach_yellow'
    'pear'
    'pearl'
    'pearl_aqua'
    'pearly_purple'
    'peridot'
    'periwinkle'
    'persian_blue'
    'persian_green'
    'persian_indigo'
    'persian_orange'
    'persian_pink'
    'persian_plum'
    'persian_red'
    'persian_rose'
    'persimmon'
    'peru'
    'phlox'
    'phthalo_blue'
    'phthalo_green'
    'piggy_pink'
    'pine_green'
    'pink'
    'pink_lace'
    'pink_orange'
    'pink_pearl'
    'pink_sherbet'
    'pistachio'
    'platinum'
    'plum_traditional'
    'plum_web'
    'portland_orange'
    'powder_blue_web'
    'princeton_orange'
    'prune'
    'prussian_blue'
    'psychedelic_purple'
    'puce'
    'pumpkin'
    'purple_heart'
    'purple_html_css'
    'purple_mountain_majesty'
    'purple_munsell'
    'purple_pizzazz'
    'purple_taupe'
    'purple_x11'
    'quartz'
    'rackley'
    'radical_red'
    'rajah'
    'raspberry'
    'raspberry_glace'
    'raspberry_pink'
    'raspberry_rose'
    'raw_umber'
    'razzle_dazzle_rose'
    'razzmatazz'
    'red'
    'red_brown'
    'red_devil'
    'red_munsell'
    'red_ncs'
    'red_orange'
    'red_pigment'
    'red_ryb'
    'red_violet'
    'redwood'
    'regalia'
    'resolution_blue'
    'rich_black'
    'rich_brilliant_lavender'
    'rich_carmine'
    'rich_electric_blue'
    'rich_lavender'
    'rich_lilac'
    'rich_maroon'
    'rifle_green'
    'robin_egg_blue'
    'rose'
    'rose_bonbon'
    'rose_ebony'
    'rose_gold'
    'rose_madder'
    'rose_pink'
    'rose_quartz'
    'rose_taupe'
    'rose_vale'
    'rosewood'
    'rosso_corsa'
    'rosy_brown'
    'royal_azure'
    'royal_blue_traditional'
    'royal_blue_web'
    'royal_fuchsia'
    'royal_purple'
    'royal_yellow'
    'rubine_red'
    'ruby'
    'ruby_red'
    'ruddy'
    'ruddy_brown'
    'ruddy_pink'
    'rufous'
    'russet'
    'rust'
    'rusty_red'
    'sacramento_state_green'
    'saddle_brown'
    'safety_orange_blaze_orange'
    'saffron'
    'salmon'
    'salmon_pink'
    'sand'
    'sand_dune'
    'sandstorm'
    'sandy_brown'
    'sandy_taupe'
    'sangria'
    'sap_green'
    'sapphire'
    'sapphire_blue'
    'satin_sheen_gold'
    'scarlet'
    'scarlet_crayola'
    'school_bus_yellow'
    'screamin_green'
    'sea_blue'
    'sea_green'
    'seal_brown'
    'seashell'
    'selective_yellow'
    'sepia'
    'shadow'
    'shamrock_green'
    'shocking_pink'
    'shocking_pink_crayola'
    'sienna'
    'silver'
    'sinopia'
    'skobeloff'
    'sky_blue'
    'sky_magenta'
    'slate_blue'
    'slate_gray'
    'smalt_dark_powder_blue'
    'smokey_topaz'
    'smoky_black'
    'snow'
    'spiro_disco_ball'
    'spring_bud'
    'spring_green'
    'st_patrick_s_blue'
    'steel_blue'
    'stil_de_grain_yellow'
    'stizza'
    'stormcloud'
    'straw'
    'sunglow'
    'sunset'
    'tan'
    'tangelo'
    'tangerine'
    'tangerine_yellow'
    'tango_pink'
    'taupe'
    'taupe_gray'
    'tea_green'
    'tea_rose_orange'
    'tea_rose_rose'
    'teal'
    'teal_blue'
    'teal_green'
    'telemagenta'
    'tenn_tawny'
    'terra_cotta'
    'thistle'
    'thulian_pink'
    'tickle_me_pink'
    'tiffany_blue'
    'tiger_s_eye'
    'timberwolf'
    'titanium_yellow'
    'tomato'
    'toolbox'
    'topaz'
    'tractor_red'
    'trolley_grey'
    'tropical_rain_forest'
    'true_blue'
    'tufts_blue'
    'tumbleweed'
    'turkish_rose'
    'turquoise'
    'turquoise_blue'
    'turquoise_green'
    'tuscan_red'
    'twilight_lavender'
    'tyrian_purple'
    'ua_blue'
    'ua_red'
    'ube'
    'ucla_blue'
    'ucla_gold'
    'ufo_green'
    'ultra_pink'
    'ultramarine'
    'ultramarine_blue'
    'umber'
    'unbleached_silk'
    'united_nations_blue'
    'university_of_california_gold'
    'unmellow_yellow'
    'up_forest_green'
    'up_maroon'
    'upsdell_red'
    'urobilin'
    'usafa_blue'
    'usc_cardinal'
    'usc_gold'
    'utah_crimson'
    'vanilla'
    'vegas_gold'
    'venetian_red'
    'verdigris'
    'vermilion_cinnabar'
    'vermilion_plochere'
    'veronica'
    'violet'
    'violet_blue'
    'violet_color_wheel'
    'violet_ryb'
    'violet_web'
    'viridian'
    'vivid_auburn'
    'vivid_burgundy'
    'vivid_cerise'
    'vivid_tangerine'
    'vivid_violet'
    'warm_black'
    'waterspout'
    'wenge'
    'wheat'
    'white'
    'white_smoke'
    'wild_blue_yonder'
    'wild_strawberry'
    'wild_watermelon'
    'wine'
    'wine_dregs'
    'wisteria'
    'wood_brown'
    'xanadu'
    'yale_blue'
    'yellow'
    'yellow_green'
    'yellow_munsell'
    'yellow_ncs'
    'yellow_orange'
    'yellow_process'
    'yellow_ryb'
    'zaffre'
    'zinnwaldite_brown'};

RGBs = ...
      [93         138         168   
        0          48         143   
      114         160         193   
      163          38          56   
      240         248         255   
      227          38          54   
      196          98          16   
      239         222         205   
      229          43          80   
      255         191           0   
      255         126           0   
      255           3          62   
      153         102         204   
      164         198          57   
      242         243         244   
      205         149         117   
      145          92         131   
      132          27          45   
      250         235         215   
        0         128           0   
      141         182           0   
      251         206         177   
        0         255         255   
      127         255         212   
       75          83          32   
       59          68          75   
      233         214         107   
      178         190         181   
      135         169         107   
      255         153         102   
      165          42          42   
      253         238           0   
      110         127         128   
       86         130           3   
        0         127         255   
      240         255         255   
      137         207         240   
      161         202         241   
      244         194         194   
       33         171         205   
      250         231         181   
      255         225          53   
      124          10           2   
      132         132         130   
      152         119         123   
      188         212         230   
      159         129         112   
      245         245         220   
      156          37          66   
      255         228         196   
       61          43          31   
      254         111          94   
      191          79          81   
        0           0           0   
       61          12           2   
       37          53          41   
       59          60          54   
      255         235         205   
      165         113         100   
       49         140         231   
      172         229         238   
      250         240         190   
        0           0         255   
      162         162         208   
       31         117         254   
      102         153         204   
       13         152         186   
        0         147         175   
        0         135         189   
       51          51         153   
        2          71         254   
       18          97         128   
      138          43         226   
      222          93         131   
      121          68          59   
        0         149         182   
      227         218         201   
      204           0           0   
        0         106          78   
      135          50          96   
        0         112         255   
      181         166          66   
      203          65          84   
       29         172         214   
      102         255           0   
      191         148         228   
      195          33          72   
      255           0         127   
        8         232         222   
      209         159         232   
      244         187         255   
      255          85         163   
      251          96         127   
        0          66          37   
      205         127          50   
      150          75           0   
      165          42          42   
      255         193         204   
      231         254         255   
      240         220         130   
       72           6           7   
      128           0          32   
      222         184         135   
      204          85           0   
      233         116          81   
      138          51          36   
      189          51         164   
      112          41          99   
       83         104         114   
       95         158         160   
      145         163         176   
        0         107          60   
      237         135          45   
      227           0          34   
      255         246           0   
      166         123          91   
       75          54          33   
       30          77          43   
      163         193         173   
      193         154         107   
      239         187         204   
      120         134         107   
      255         239           0   
      255           8           0   
      228         113         122   
        0         191         255   
       89          39          32   
      196          30          58   
        0         204         153   
      150           0          24   
      215           0          64   
      235          76          66   
      255           0          56   
      255         166         201   
      179          27          27   
      153         186         221   
      237         145          33   
        6          42         120   
      146         161         207   
      172         225         175   
        0         123         167   
       47         132         124   
      178         255         255   
       73         151         208   
      222          49          99   
      236          59         131   
        0         123         167   
       42          82         190   
      109         155         195   
        0         122         165   
      224          60          49   
      160         120          90   
      250         214         165   
       54          69          79   
      230         143         172   
      223         255           0   
      127         255           0   
      222          49          99   
      255         183         197   
      205          92          92   
      222         111         161   
      168          81         110   
      170          56          30   
      123          63           0   
      210         105          30   
      255         167           0   
      152         129         123   
      227          66          52   
      210         105          30   
      228         208          10   
      251         204         231   
        0          71         171   
      210         105          30   
      111          78          55   
      155         221         255   
      248         131         121   
        0          46          99   
      140         146         172   
      184         115          51   
      218         138         103   
      173         111         105   
      203         109          81   
      153         102         102   
      255          56           0   
      255         127          80   
      248         131         121   
      255          64          64   
      137          63          69   
      251         236          93   
      179          27          27   
      100         149         237   
      255         248         220   
      255         248         231   
      255         188         217   
      255         253         208   
      220          20          60   
      190           0          50   
        0         255         255   
        0         183         235   
      255         255          49   
      240         225          48   
        0           0         139   
      101          67          33   
       93          57          84   
      164           0           0   
        8          69         126   
      152         105          96   
      205          91          69   
        0         139         139   
       83         104         120   
      184         134          11   
      169         169         169   
        1          50          32   
        0          65         106   
       26          36          33   
      189         183         107   
       72          60          50   
      115          79         150   
      139           0         139   
        0          51         102   
       85         107          47   
      255         140           0   
      153          50         204   
      119         158         203   
        3         192          60   
      150         111         214   
      194          59          34   
      231          84         128   
        0          51         153   
      135          38          87   
      139           0           0   
      233         150         122   
       86           3          25   
      143         188         143   
       60          20          20   
       72          61         139   
       47          79          79   
       23         114          69   
      145         129          81   
      255         168          18   
       72          60          50   
      204          78          92   
        0         206         209   
      148           0         211   
      155         135          12   
        0         112          60   
       85          85          85   
      215          10          83   
      169          32          62   
      239          48          56   
      233         105          44   
      218          50         135   
      250         214         165   
      185          78          72   
      112          66          65   
      193          84         193   
        0          75          73   
      153          85         187   
      204           0         204   
      255         203         164   
      255          20         147   
      132          63          91   
      255         153          51   
        0         191         255   
      102          66          77   
       21          96         189   
      193         154         107   
      237         201         175   
      105         105         105   
       30         144         255   
      215          24         104   
      133         187         101   
      150         113          23   
        0           0         156   
      225         169          95   
       85          93          80   
      194         178         128   
       97          64          81   
      240         234         214   
       16          52         166   
      125         249         255   
      255           0          63   
        0         255         255   
        0         255           0   
      111           0         255   
      244         187         255   
      204         255           0   
      191           0         255   
       63           0         255   
      143           0         255   
      255         255           0   
       80         200         120   
      180         131         149   
      150         200         162   
      193         154         107   
      128          24          24   
      181          51         137   
      244           0         161   
      229         170         112   
       77          93          83   
       79         121          66   
      255          40           0   
      108          84          30   
      206          32          41   
      178          34          34   
      226          88          34   
      252         142         172   
      247         233         142   
      238         220         130   
      255         250         240   
      255         191           0   
      255          20         147   
      204         255           0   
      255           0          79   
        1          68          33   
       34         139          34   
      166         123          91   
        0         114         187   
      134          96         142   
      204         255           0   
      199          44          72   
      246          74         138   
      255           0         255   
      193          84         193   
      255         119         255   
      199          67         117   
      228         132           0   
      204         102         102   
      220         220         220   
      228         155          15   
      248         248         255   
      176         101           0   
       96         130         182   
      230         232         250   
      212         175          55   
      255         215           0   
      153         101          21   
      252         194           0   
      255         223           0   
      218         165          32   
      168         228         160   
      128         128         128   
       70          89          69   
      128         128         128   
      190         190         190   
        0         255           0   
       28         172         120   
        0         128           0   
        0         168         119   
        0         159         107   
        0         165          80   
      102         176          50   
      173         255          47   
      169         154         134   
        0         255         127   
      102          56          84   
       68         108         207   
       82          24         250   
      233         214         107   
       63         255           0   
      201           0          22   
      218         145           0   
      128         128           0   
      223         115         255   
      244           0         161   
      240         255         240   
        0         127         191   
       73         121         107   
      255          29         206   
      255         105         180   
       53          94          59   
      113         166         210   
      252         247          94   
        0          35         149   
      178         236          93   
       19         136           8   
      205          92          92   
      227         168          87   
      111           0         255   
        0          65         106   
       75           0         130   
        0          47         167   
      255          79           0   
      186          22          12   
      192          54          44   
       90          79         207   
      244         240         236   
        0         144           0   
      255         255         240   
        0         168         107   
      248         222         126   
      215          59          62   
      165          11          94   
       52          52          52   
      250         218          94   
      189         218          87   
       41         171         135   
       76         187          23   
      124          28           5   
      195         176         145   
      240         230         140   
      232           0          13   
        8         120          48   
      214         202         221   
       38          97         156   
      254         254          34   
      169         186         157   
      207          16          32   
      204         204         255   
      255         240         245   
      181         126         220   
      196         195         208   
      148          87         235   
      238         130         238   
      230         230         250   
      251         174         210   
      150         123         182   
      251         160         227   
      230         230         250   
      124         252           0   
      255         247           0   
      255         250         205   
      227         255           0   
       26          17          16   
      253         213         177   
      173         216         230   
      181         101          29   
      230         103         113   
      240         128         128   
      147         204         234   
      245         105         145   
      224         255         255   
      249         132         239   
      250         250         210   
      211         211         211   
      144         238         144   
      240         230         140   
      177         156         217   
      255         182         193   
      233         116          81   
      255         160         122   
      255         153         153   
       32         178         170   
      135         206         250   
      119         136         153   
      179         139         109   
      230         143         172   
      255         255         224   
      200         162         200   
      191         255           0   
       50         205          50   
        0         255           0   
      157         194           9   
       25          89           5   
      250         240         230   
      193         154         107   
      108         160         220   
       83          75          79   
      230          32          32   
      255           0         255   
      202          31         123   
      255           0         144   
      170         240         209   
      248         244         255   
      192          64           0   
      251         236          93   
       96          80         220   
       11         218          81   
      151         154         170   
      255         130          67   
      116         195         101   
      136           0         133   
      195          33          72   
      128           0           0   
      176          48          96   
      224         176         255   
      145          95         109   
      239         152         170   
      115         194         251   
      229         183          59   
      102         221         170   
        0           0         205   
      226           6          44   
      175          64          53   
      243         229         171   
        3          80         150   
       28          53          45   
      221         160         221   
      186          85         211   
        0         103         165   
      147         112         219   
      187          51         133   
      170          64         105   
       60         179         113   
      123         104         238   
      201         220         135   
        0         250         154   
      103          76          71   
       72         209         204   
      121          68          59   
      217          96          59   
      199          21         133   
      248         184         120   
      248         222         126   
      253         188         180   
       25          25         112   
        0          73          83   
      255         196          12   
       62         180         137   
      245         255         250   
      152         255         152   
      255         228         225   
      250         235         215   
      150         113          23   
      115         169         194   
      174          12           0   
      173         223         173   
       48         186         143   
      153         122         141   
       24          69          59   
      197          75         140   
      255         219          88   
       33          66          30   
      246         173         198   
       42         128           0   
      250         218          94   
      255         222         173   
        0           0         128   
      255         163          67   
      254          65         100   
       57         255          20   
      215         131         127   
      164         221         237   
        5         144          51   
        0         119         190   
      204         119          34   
        0         128           0   
      207         181          59   
      253         245         230   
      121         104         120   
      103          49          71   
      192         128         129   
      128         128           0   
       60          52          31   
      107         142          35   
      154         185         115   
       53          56          57   
      183         132         167   
      255         127           0   
      255         159           0   
      255          69           0   
      251         153           2   
      255         165           0   
      218         112         214   
      101          67          33   
      153           0           0   
       65          74          76   
      255         110          74   
        0          33          71   
        0         102           0   
       39          59         226   
      104          40          96   
      188         212         230   
      175         238         238   
      152         118          84   
      175          64          53   
      155         196         226   
      221         173         175   
      218         138         103   
      171         205         239   
      230         190         138   
      238         232         170   
      152         251         152   
      220         208         255   
      249         132         229   
      250         218         221   
      221         160         221   
      219         112         147   
      150         222         209   
      201         192         187   
      236         235         189   
      188         152         126   
      219         112         147   
      120          24          74   
      255         239         213   
       80         200         120   
      174         198         207   
      131         105          83   
      207         207         196   
      119         221         119   
      244         154         194   
      255         179          71   
      222         165         164   
      179         158         181   
      255         105          97   
      203         153         201   
      253         253         150   
      128           0         128   
       83         104         120   
      255         229         180   
      255         203         164   
      255         204         153   
      255         218         185   
      250         223         173   
      209         226          49   
      234         224         200   
      136         216         192   
      183         104         162   
      230         226           0   
      204         204         255   
       28          57         187   
        0         166         147   
       50          18         122   
      217         144          88   
      247         127         190   
      112          28          28   
      204          51          51   
      254          40         162   
      236          88           0   
      205         133          63   
      223           0         255   
        0          15         137   
       18          53          36   
      253         221         230   
        1         121         111   
      255         192         203   
      255         221         244   
      255         153         102   
      231         172         207   
      247         143         167   
      147         197         114   
      229         228         226   
      142          69         133   
      221         160         221   
      255          90          54   
      176         224         230   
      255         143           0   
      112          28          28   
        0          49          83   
      223           0         255   
      204         136         153   
      255         117          24   
      105          53         156   
      128           0         128   
      150         120         182   
      159           0         197   
      254          78         218   
       80          64          77   
      160          32         240   
       81          72          79   
       93         138         168   
      255          53          94   
      251         171          96   
      227          11          93   
      145          95         109   
      226          80         152   
      179          68         108   
      130         102          68   
      255          51         204   
      227          37         107   
      255           0           0   
      165          42          42   
      134           1          17   
      242           0          60   
      196           2          51   
      255          83          73   
      237          28          36   
      254          39          18   
      199          21         133   
      171          78          82   
       82          45         128   
        0          35         135   
        0          64          64   
      241         167         254   
      215           0          64   
        8         146         208   
      167         107         207   
      182         102         210   
      176          48          96   
       65          72          51   
        0         204         204   
      255           0         127   
      249          66         158   
      103          72          70   
      183         110         121   
      227          38          54   
      255         102         204   
      170         152         169   
      144          93          93   
      171          78          82   
      101           0          11   
      212           0           0   
      188         143         143   
        0          56         168   
        0          35         102   
       65         105         225   
      202          44         146   
      120          81         169   
      250         218          94   
      209           0          86   
      224          17          95   
      155          17          30   
      255           0          40   
      187         101          40   
      225         142         150   
      168          28           7   
      128          70          27   
      183          65          14   
      218          44          67   
        0          86          63   
      139          69          19   
      255         103           0   
      244         196          48   
      255         140         105   
      255         145         164   
      194         178         128   
      150         113          23   
      236         213          64   
      244         164          96   
      150         113          23   
      146           0          10   
       80         125          42   
       15          82         186   
        0         103         165   
      203         161          53   
      255          36           0   
      253          14          53   
      255         216           0   
      118         255         122   
        0         105         148   
       46         139          87   
       50          20          20   
      255         245         238   
      255         186           0   
      112          66          20   
      138         121          93   
        0         158          96   
      252          15         192   
      255         111         255   
      136          45          23   
      192         192         192   
      203          65          11   
        0         116         116   
      135         206         235   
      207         113         175   
      106          90         205   
      112         128         144   
        0          51         153   
      147          61          65   
       16          12           8   
      255         250         250   
       15         192         252   
      167         252           0   
        0         255         127   
       35          41         122   
       70         130         180   
      250         218          94   
      153           0           0   
       79         102         106   
      228         217         111   
      255         204          51   
      250         214         165   
      210         180         140   
      249          77           0   
      242         133           0   
      255         204           0   
      228         113         122   
       72          60          50   
      139         133         137   
      208         240         192   
      248         131         121   
      244         194         194   
        0         128         128   
       54         117         136   
        0         130         127   
      207          52         118   
      205          87           0   
      226         114          91   
      216         191         216   
      222         111         161   
      252         137         172   
       10         186         181   
      224         141          60   
      219         215         210   
      238         230           0   
      255          99          71   
      116         108         192   
      255         200         124   
      253          14          53   
      128         128         128   
        0         117          94   
        0         115         207   
       65         125         193   
      222         170         136   
      181         114         129   
       48         213         200   
        0         255         239   
      160         214         180   
      124          72          72   
      138          73         107   
      102           2          60   
        0          51         170   
      217           0          76   
      136         120         195   
       83         104         149   
      255         179           0   
       60         208         112   
      255         111         255   
       18          10         143   
       65         102         245   
       99          81          71   
      255         221         202   
       91         146         229   
      183         135          39   
      255         255         102   
        1          68          33   
      123          17          19   
      174          32          41   
      225         173          33   
        0          79         152   
      153           0           0   
      255         204           0   
      211           0          63   
      243         229         171   
      197         179          88   
      200           8          21   
       67         179         174   
      227          66          52   
      217          96          59   
      160          32         240   
      143           0         255   
       50          74         178   
      127           0         255   
      134           1         175   
      238         130         238   
       64         130         109   
      146          39          36   
      159          29          53   
      218          29         129   
      255         160         137   
      159           0         255   
        0          66          66   
      164         244         249   
      100          84          82   
      245         222         179   
      255         255         255   
      245         245         245   
      162         173         208   
      255          67         164   
      252         108         133   
      114          47          55   
      103          49          71   
      201         160         220   
      193         154         107   
      115         134         120   
       15          77         146   
      255         255           0   
      154         205          50   
      239         204           0   
      255         211           0   
      255         174          66   
      255         239           0   
      254         254          51   
        0          20         168   
       44          22           8];

RGBs = num2cell(RGBs,2);

% Make a hash map for O(1) search
map = containers.Map(colors,RGBs);
if map.isKey(color)
    rgb = map(color) ./ 255;
else
    rgb = [0 0 0];
end

end
