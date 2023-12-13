#' @export

getEOM <-function(){
  eom <- pops.sdm::state(c('Wisconsin', 'Michigan',
                           'Illinois', 'Indiana', 'Ohio',
                           'Kentucky', 'Tennessee',
                           'Mississippi', 'Alabama',
                           'Florida', 'Georgia',
                           'South Carolina', 'North Carolina',
                           'Virginia', 'West Virgina',
                           'Maryland', 'District of Columbia',
                           'Delaware', 'Pennsylvania',
                           'New Jersey', 'New York',
                           'Connecticut', 'Rhode Island',
                           'Massachusetts', 'Vermont',
                           'New Hampshire', 'Maine'
  ))
  return(eom)
}

getEOR <-function(){
  eor <- pops.sdm::state(c('North Dakota', 'South Dakota',
                           'Nebraska', 'Kansas', 'Arkansas',
                           'Oklahoma', 'Texas', 'Louisiana',
                           'Missouri', 'Iowa',
                           'Minnesota', 'Wisconsin', 'Michigan',
                           'Illinois', 'Indiana', 'Ohio',
                           'Kentucky', 'Tennessee',
                           'Mississippi', 'Alabama',
                           'Florida', 'Georgia',
                           'South Carolina', 'North Carolina',
                           'Virginia', 'West Virginia',
                           'Maryland', 'District of Columbia',
                           'Delaware', 'Pennsylvania',
                           'New Jersey', 'New York',
                           'Connecticut', 'Rhode Island',
                           'Massachusetts', 'Vermont',
                           'New Hampshire', 'Maine'
  ))
  return(eor)
}

getWOM <- function(){
  wom <- pops.sdm::state(c('Washington', 'Oregon', 'California',
                           'Arizona', 'Nevada', 'Utah', 'Idaho',
                           'Montana', 'Wyoming', 'Colorado',
                           'New Mexico', 'Texas', 'Oklahoma',
                           'Kansas', 'Nebraska', 'South Dakota',
                           'North Dakota', 'Minnesota', 'Iowa',
                           'Missouri', 'Arkansas', 'Loiusiana'
  ))
  return(wom)
}

getWOR <- function(){
  wor <- pops.sdm::state(c('Washington', 'Oregon', 'California',
                           'Arizona', 'Nevada', 'Utah', 'Idaho',
                           'Montana', 'Wyoming', 'Colorado', 'New Mexico'
  ))
  return(wor)
}

getMILP <- function(){
  milp <- pops.sdm::county(state='Michigan', names=c('Alcona', 'Montmorency', 'Oscoda', 'Cass',
                                                     'Iosco', 'Antrim', 'Wayne', 'Branch',
                                                     'Alpena', 'Arenac', 'Benzie', 'Crawford',
                                                     'Grand Traverse', 'Gratiot', 'Ingham',
                                                     'Eaton', 'Charlevoix', 'Barry', 'Kent',
                                                     'Emmet', 'Genesee', 'Berrien', 'Huron',
                                                     'Calhoun', 'Hillsdale', 'Ionia', 'Monroe',
                                                     'Kalkaska', 'Cheboygan', 'Wexford', 'Bay',
                                                     'Missaukee', 'Presque Isle', 'Washtenaw',
                                                     'Shiawassee', 'Oakland', 'Livingston',
                                                     'Tuscola', 'Montcalm', 'Lenawee', 'Ogemaw',
                                                     'St. Clair',  'Macomb', 'Sanilac', 'Mason',
                                                     'Otsego', 'Manistee', 'Roscommon', 'Lake',
                                                     'Osceola', 'Lapeer', 'Gladwin', 'Oceana',
                                                     'Newaygo', 'Mecosta', 'Midland', 'Ottawa',
                                                     'Clinton', 'Muskegon', 'Saginaw', 'Clare',
                                                     'Allegan', 'Van Buren', 'St. Joseph',
                                                     'Kalamazoo', 'Jackson', 'Isabella'))
  return(milp)
}
