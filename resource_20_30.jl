# Parameters of the resource allocation model with 20 resources and 30 customer types (adapted from Luedtke, 2014)

const numResources = 20 
const numCustomers = 30 

# first-stage objective coefficients
const costVector = [0.7432683327183562 1.1091782740319707 1.153303321756463 0.7751143229777455 0.8098221660821516 1.1324011945408334 0.9859746207969543 0.9502566702752542 1.1909969906693632 0.9630816483432164 1.200255951169507 1.1412400075515754 1.0668906064783976 0.6939486330511746 1.3054592887126564 1.121268771178427 1.0238773621726285 0.9083077835585992 1.2091835133144502 1.171512305637218 ] 

# recourse objective coefficients
const recourseCostCoeff = [2.3130311834976878 1.970110451642846 2.079834948369165 2.180997535033276 2.1887583154880983 2.157134900936415 2.148981366909944 2.063591377530072 2.114553485889044 2.2016770965819905 1.9966732632984514 2.2330830064352423 2.073101936010451 2.078736481418292 2.1023788717803735 2.20322688021812 2.0598950590963474 2.186612946072563 2.2213945559591637 2.4194867761920977 2.177263511555545 2.2031710809006455 1.9909192466456946 2.2396747891757802 2.1661606250152716 2.049535136705102 2.1171635559467483 2.2147607076134714 2.203825916684725 2.3108454507501515 ] 

# yield parameters
const rho = [0.9555664467976429 0.9739615542635913 0.9122667177690378 0.9340610389388502 0.9246641057347171 0.9945151890629177 0.9615143263620417 0.9659480369131807 0.9358045688085963 0.9147377425339717 0.9756453545272545 0.9350544644953063 0.9839239868461557 0.9953649754312432 0.9288115703035881 0.9946050201412094 0.9575793792880826 0.9224016204902284 0.905702317419421 0.9507354965994901 ] 

# service rate parameters
const mu = zeros(Float64,20,30) 
mu[1,:] = [0.0 1.5425864608494138 1.8677174853336527 1.7036521533592883 0.0 0.0 1.7940537688365086 1.8901902999735647 0.0 1.743223975958845 0.0 0.0 1.9010265875424377 1.8193631099999497 0.0 0.0 0.0 0.0 1.7774304885999672 0.0 1.9956705204853944 1.7673613846413179 0.0 1.6465508373626938 1.7811853533294304 0.0 1.7386645478259126 0.0 0.0 1.5126998082534597 ] 
mu[2,:] = [2.2607279588202296 1.9084964021630282 0.0 2.0695620946729028 2.4463795504578956 0.0 2.159963710150123 2.256100241287179 0.0 2.1091339172724597 2.2022354542170515 0.0 0.0 2.185273051313564 1.990070719163648 0.0 0.0 2.27455757738075 2.1433404299135814 0.0 2.361580461799009 2.133271325954932 2.2274657922865684 0.0 2.147095294643045 2.3868225163693566 0.0 2.0854662478744874 2.197861409557054 0.0 ] 
mu[3,:] = [0.0 1.9526214498875207 0.0 0.0 2.490504598182388 2.12721997941551 2.204088757874615 2.3002252890116717 2.362568160228464 2.153258964996952 2.2463605019415436 0.0 2.3110615765805447 0.0 2.03419576688814 0.0 2.073572838292769 0.0 2.187465477638074 1.8667986470938884 2.405705509523501 2.1773963736794246 2.2715908400110605 0.0 2.1912203423675374 2.430947564093849 0.0 2.1295912955989795 0.0 0.0 ] 
mu[4,:] = [1.9266640077660047 1.574432451108803 1.899563475593042 1.7354981436186776 2.1123155994036704 0.0 1.8258997590958979 1.922036290232954 0.0 1.7750699662182343 1.8681715031628263 0.0 0.0 1.851209100259339 1.6560067681094228 0.0 0.0 0.0 0.0 1.488609648315171 2.0275165107447837 1.7992073749007071 0.0 0.0 0.0 2.0527585653151315 0.0 1.751402296820262 1.863797458502829 1.544545798512849 ] 
mu[5,:] = [0.0 0.0 0.0 0.0 2.1470234425080768 1.7837388237411989 1.860607602200304 1.95674413333736 0.0 1.8097778093226404 0.0 0.0 1.967580420906233 1.885916943363745 0.0 1.8935808138490131 0.0 0.0 1.8439843219637626 1.523317491419577 2.0622243538491896 0.0 1.9281096843367493 0.0 0.0 2.087466408419538 0.0 1.7861101399246682 1.898505301607235 0.0 ] 
mu[6,:] = [2.2839508793290926 0.0 2.2568503471561296 2.0927850151817653 2.4696024709667586 0.0 2.1831866306589855 0.0 2.3416660330128343 2.1323568377813222 0.0 0.0 2.290159449364915 2.208495971822427 2.0132936396725105 0.0 0.0 2.2977804978896126 2.1665633504224444 0.0 2.3848033823078714 0.0 2.250688712795431 0.0 2.170318215151908 2.4100454368782196 0.0 2.10868916838335 0.0 1.9018326700759371 ] 
mu[7,:] = [2.1375243055852136 0.0 2.1104237734122506 0.0 2.3231758972228795 0.0 2.0367600569151065 2.1328965880521626 0.0 1.9859302640374432 2.079031800982035 1.806567963946108 2.1437328756210356 0.0 1.8668670659286315 0.0 1.9062441373332604 0.0 2.0201367766785654 1.6994699461343798 0.0 2.010067672719916 0.0 0.0 0.0 2.2636188631343406 0.0 1.962262594639471 2.0746577563220376 0.0 ] 
mu[8,:] = [0.0 1.7495747984063117 2.0747058228905506 1.9106404909161863 0.0 0.0 2.0010421063934065 2.0971786375304626 0.0 0.0 0.0 0.0 0.0 2.0263514475568476 1.8311491154069315 2.0340153180421154 0.0 0.0 0.0 0.0 2.2026588580422923 1.9743497221982158 0.0 1.8535391749195917 1.9881736908863283 2.2279009126126406 0.0 0.0 2.0389398058003376 1.7196881458103577 ] 
mu[9,:] = [2.3425466754576223 1.990315118800421 2.3154461432846594 2.151380811310295 0.0 2.1649136483284104 2.2417824267875153 2.337918957924572 2.400261829141364 2.190952633909852 0.0 2.011590333818517 2.348755245493445 2.267091767950957 2.0718894358010402 2.2747556384362246 2.111266507205669 2.3563762940181423 0.0 0.0 2.443399178436401 2.215090042592325 2.3092845089239606 0.0 2.2289140112804375 0.0 2.1863932057769193 0.0 0.0 1.9604284662044669 ] 
mu[10,:] = [0.0 1.762399776474274 2.087530800958513 1.9234654689841484 2.3002829247691414 1.9369983060022635 2.013867084461369 0.0 2.1723464868152176 1.963037291583705 2.0561388285282973 1.78367499149237 2.120839903167298 2.03917642562481 0.0 0.0 0.0 0.0 1.9972438042248273 1.6765769736806417 2.2154838361102547 0.0 2.0813691665978142 0.0 2.0009986689542907 2.2407258906806025 0.0 0.0 2.0517647838683 0.0 ] 
mu[11,:] = [2.351805635957766 0.0 0.0 2.160639771810439 0.0 0.0 2.2510413872876596 2.3471779184247152 2.4095207896415083 0.0 2.293313131354588 0.0 0.0 0.0 2.0811483963011845 0.0 0.0 0.0 0.0 0.0 2.4526581389365454 2.2243490030924686 2.318543469424105 2.103538455813845 2.238172971780581 2.477900193506893 2.195652166277063 2.1765439250120235 2.2889390866945902 0.0 ] 
mu[12,:] = [2.2927896923398343 0.0 0.0 2.1016238281925075 2.4784412839775003 0.0 2.1920254436697277 0.0 0.0 2.1411956507920644 2.234297187736656 0.0 2.298998262375657 0.0 2.0221324526832527 2.224998655318437 2.0615095240878816 2.3066193109003548 0.0 1.8547353328890008 2.3936421953186136 0.0 2.259527525806173 0.0 2.1791570281626496 2.4188842498889613 0.0 2.117527981394092 2.229923143076659 1.910671483086679 ] 
mu[13,:] = [2.2184402912666568 0.0 0.0 2.0272744271193295 0.0 0.0 2.1176760425965497 2.2138125737336063 0.0 0.0 0.0 1.8874839496275513 2.2246488613024793 2.1429853837599913 1.947783051610075 2.150649254245259 1.9871601230147036 2.2322699098271768 2.1010527623600086 1.780385931815823 0.0 2.0909836584013592 2.185178124732995 0.0 0.0 2.344534848815784 0.0 2.043178580320914 0.0 1.8363220820135013 ] 
mu[14,:] = [0.0 1.4932667611822321 1.818397785666471 1.6543324536921067 0.0 0.0 1.744734069169327 0.0 1.9032134715231757 1.6939042762916634 0.0 1.5145419762003283 1.851706887875256 1.770043410332768 1.574841078182852 0.0 0.0 1.8593279363999538 0.0 1.4074439583886 1.9463508208182128 1.7180416849741362 0.0 0.0 1.7318656536622488 1.9715928753885608 1.689344848158731 1.6702366068936911 1.782631768576258 1.463380108586278 ] 
mu[15,:] = [2.4570089735009155 0.0 0.0 2.2658431093535882 0.0 2.2793759463717036 2.356244724830809 2.4523812559678646 2.5147241271846577 0.0 0.0 2.12605263186181 0.0 0.0 0.0 2.389217936479518 2.2257288052489623 2.4708385920614355 0.0 2.018954614050082 0.0 2.329552340635618 2.4237468069672543 2.208741793356994 2.3433763093237303 0.0 0.0 2.281747262555173 2.3941424242377396 2.07489076424776 ] 
mu[16,:] = [0.0 1.9205868993094848 2.2457179237937233 2.081652591819359 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.9418621143275807 2.2790270260025087 2.1973635484600207 0.0 0.0 0.0 0.0 2.155430927060038 1.8347640965158525 2.373670958945465 2.1453618231013887 2.2395562894330245 0.0 2.1591857917895014 2.398913013515813 0.0 0.0 0.0 0.0 ] 
mu[17,:] = [2.1754270469608876 1.8231954903036862 2.148326514787925 1.9842611828135603 0.0 1.9977940198316757 2.074662798290781 2.1707993294278367 2.2331422006446298 2.0238330054131173 0.0 1.8444707053217821 2.1816356169967097 2.0999721394542217 1.9047698073043058 2.10763600993949 0.0 2.1892566655214076 2.0580395180542395 0.0 0.0 0.0 0.0 1.9271598668169663 2.0617943827837024 2.3015216045100146 0.0 0.0 2.1125604976977117 1.7933088377077322 ] 
mu[18,:] = [2.0598574683468582 1.7076259116896568 0.0 0.0 2.245509059984524 0.0 1.9590932196767517 2.0552297508138078 2.1175726220306004 1.908263426799088 2.00136496374368 1.7289011267077528 2.0660660383826808 1.9844025608401927 1.7892002286902766 1.9920664313254606 1.828577300094905 0.0 1.94246993944021 1.6218031088960245 2.1607099713256375 1.9324008354815607 2.026595301813197 0.0 0.0 0.0 0.0 1.8845957574011156 1.9969909190836828 0.0 ] 
mu[19,:] = [2.360733198102709 2.0085016414455077 0.0 2.1695673339553823 0.0 2.183100170973497 0.0 0.0 0.0 0.0 0.0 2.029776856463604 2.3669417681385316 2.2852782905960436 2.0900759584461275 2.292942161081312 2.1294530298507564 2.3745628166632295 0.0 1.9226788386518756 2.4615857010814883 2.2332765652374116 2.327471031569048 2.1124660179587877 2.2471005339255243 0.0 2.2045797284220066 2.185471487156967 2.2978666488395336 0.0 ] 
mu[20,:] = [0.0 1.9708304337682758 2.2959614582525143 2.13189612627815 2.508713582063143 0.0 2.22229774175537 2.3184342728924268 2.380777144109219 0.0 2.2645694858222987 1.9921056487863718 0.0 0.0 2.052404750768895 2.2552709534040796 2.091781822173524 2.3368916089859972 2.205674461518829 0.0 2.423914493404256 0.0 2.2897998238918156 0.0 2.2094293262482925 0.0 2.1669085207447742 2.1478002794797346 2.2601954411623018 0.0 ] 