    rm(list = ls())

    setwd("C:/Users/ezraa/OneDrive/Desktop/GPCR mutations/Alpha subunits/Binding partners")

    library(dplyr)
    library(tidyr)
    library(stringr)
    library(httr)
    library(rvest)


directory = "C:/Users/ezraa/OneDrive/Desktop/GPCR mutations/Alpha subunits/Binding partners/"
    get_bonds = function(url){

        retrieve_url = GET(url) 
        text = content(retrieve_url, as = "text") %>% str_split(pattern = "\n") 
        l = text[[1]]

    l = readLines(url)

    empty_df = data.frame(matrix(ncol = 14))

    search_l = function(start, stop, bond){
        data = l[(grep(pattern = start,l)+7):(grep(pattern=stop,l)-2)] %>%
            lapply(trimws) %>% 
            lapply(FUN = gsub, pattern = "\\s+", replacement = ",") %>% 
            mapply(FUN = str_split,",") %>% 
            as.data.frame() %>% 
            t() %>% 
            as.data.frame() %>% 
            mutate(bond = bond)
        return(data)
    }

    search_l_end = function(start,bond) {
        l[((grep(pattern = start,l))+7):(grep(l, pattern = "Number of")[1]-2)]  %>% 
        lapply(trimws) %>% 
        lapply(FUN = gsub, pattern = "\\s+", replacement = ",") %>% 
        mapply(FUN = str_split,",") %>% 
        as.data.frame() %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate(bond = bond)
    }

    if(grep(pattern = "Hydrogen bonds", l) %>% length() == 0) {
        hbs = empty_df
        } else if (grep(pattern = "Non-bonded contacts", l) %>% length() != 0) {
        hbs = search_l(start = "Hydrogen bonds", stop = "Non-bonded contacts", bond = "HBs")
        } else if (grep(pattern = "Salt bridges", l) %>% length() != 0) {
        hbs = search_l(start = "Hydrogen bonds", stop = "Salt bridges", bond = "HBs")
        } else {
        hbs = search_l_end(start = "Hydrogen bonds", bond = "HBs")
        }

    if(grep(pattern = "Non-bonded contacts", l) %>% length() == 0) {
        nbcs = empty_df 
        } else if (grep(pattern = "Salt bridges", l) %>% length() != 0) {
        nbcs = search_l(start = "Non-bonded contacts", stop = "Salt bridges", bond = "SBCs")
        } else {
        nbcs = search_l_end(start = "Non-bonded contacts", bond = "SBCs")
        }


    if(grep(pattern = "Salt bridges",l) %>% length() == 0){
        sbs = empty_df 
        } else {
            sbs = search_l_end(start = "Salt bridges", bond = "SBs")
            
        }


    # get pdb code
    pdb.line = l[grep(pattern = "PDB code", l)]
    regex = "(?<=PDB code:).*?(?=Chains)"
    matches = regmatches(pdb.line, gregexpr(regex,pdb.line, perl = TRUE)) %>% paste() %>% trimws()

    cols = c("Number", "Atom 1 no.", "Atom 1 name", "Res 1 name", "Res 1 no.", "Chain 1", "Arrow", "Atom 2 no.", "Atom 2 name", "Res 2 name", "Res 2 no.", "Chain 2", "Distance", "bond")

    colnames(hbs) = cols
    colnames(nbcs) = cols
    colnames(sbs) = cols  

    all = rbind(hbs,nbcs,sbs)  %>% mutate(pdb.code = matches)

    chain1 = str_match(string = pdb.line, pattern = "Chains(.*?)\\}")[2] %>% trimws()

    chain2 = str_match(string = pdb.line, pattern = "\\{(.*)")[,2] %>% trimws()

    # filename = paste(
    #     "generated CSVs/", matches, " ", chain1, " to ", chain2, ".csv",
    #     sep = ""
    # )

    return(all)
    }

    chain1 = "A"
    chain2 = "B"
    id = "7jsn"
    bonds = function(id, chain1,chain2, directory){
        url_text = paste("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=",id,"&chain1=",chain1,"&chain2=",chain2,
        sep = "")
        df = get_bonds(url_text)

        chain1_url = paste("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=7jsn&template=protein.html&r=wiring&l=1&chain=",chain1, sep = "")
        chain1_page = read_html(chain1_url) %>% html_nodes("tr")
        chain1_text = chain1_page %>% html_text() %>% trimws() %>% str_split("\n") %>% unlist() %>% sapply(FUN = gsub,pattern = "\\s+", replacement = "") %>% unlist()
        chain1_text = chain1_text[nzchar(chain1_text)]
        uniprot_code1 = chain1_text[grep(pattern = "UniProtcode", x = chain1_text)[1]+1]
        uniprot_id1 = chain1_text[grep(pattern = "UniProtcode", x = chain1_text)[1]+2]

        chain2_url = paste("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=7jsn&template=protein.html&r=wiring&l=1&chain=",chain2, sep = "")
        chain2_page = read_html(chain2_url) %>% html_nodes("tr")
        chain2_text = chain2_page %>% html_text() %>% trimws() %>% str_split("\n") %>% unlist() %>% sapply(FUN = gsub,pattern = "\\s+", replacement = "") %>% unlist()
        chain2_text = chain2_text[nzchar(chain2_text)]
        uniprot_code2 = chain2_text[grep(pattern = "UniProtcode", x = chain2_text)[1]+1]
        uniprot_id2 = chain2_text[grep(pattern = "UniProtcode", x = chain2_text)[1]+2]

        df = df %>% mutate(
            chain_1_id = uniprot_code1,
            chain_1_name = uniprot_id1,
            chain_2_id = uniprot_id2,
            chain_2_name = uniprot_code2,
        )


        matches = df$pdb.code[1]

        filename = paste(
            directory, matches, " ", chain1, " to ", chain2, ".csv",
            sep = ""
        )

        write.csv(df, filename)

        print(filename)
    }


    # 7f6h
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=7f6h&chain1=B&chain2=C")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=7f6h&chain1=A&chain2=B")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=7f6h&chain1=C&chain2=D")
    # bonds("7f6h","A","B")
    # bonds("7f6h","B","C")
    # bonds("76fh", "C", "D")

    # # 7sq2
    # # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=7sq2&chain1=A&chain2=B")
    # bonds("7sq2", "A", "B")

    # 60ij
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6oij&chain1=A&chain2=B")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6oij&chain1=A&chain2=H")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6oij&chain1=A&chain2=R")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6oij&chain1=B&chain2=G")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6oij&chain1=B&chain2=H")
    # get_bonds("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb=6oij&chain1=B&chain2=R")
    # bonds("60ij", "A","B")
    # bonds("60ij", "A","H")
    # bonds("60ij", "A","R")
    # bonds("60ij", "B","G")
    # bonds("60ij", "B","H")
    # bonds("60ij", "B","R")

    #7jsn
    # bonds("7jsn", "A","B")
    # bonds("7jsn", "A", "C")
    # bonds("7jsn", "A","D")
    # bonds("7jsn", "A", "F")
    # bonds("7jsn", "B","C")
    # bonds("7jsn", "B", "D")
    # bonds("7jsn", "B", "E")
    # bonds("7jsn", "C","E")
    # bonds("7jsn", "D", "F")

    # 37ck
    # https://www.rcsb.org/structure/3C7K
    # bonds("37ck", "A", "B")
    # bonds("37ck", "C", "D")

    # # 6oik 
    # bonds("60ik", "A", "B")
    # bonds("60ik", "A", "H")
    # bonds("60ik", "A", "R")
    # bonds("60ik", "B", "G")
    # bonds("60ik", "B", "H")

    # 8x2k 
    # bonds("8X2K", "A","B")
    # bonds("8X2K", "A","C")
    # bonds("8X2K", "B","C")
    # bonds("8X2K", "B","E")
    # bonds("8X2K", "C","D")
    # bonds("8X2K", "C","E")

    # 8x79
    # https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=8x79&template=interfaces.html&c=999
    # bonds("8x79","A","B")
    # bonds("8x79","A","N")
    # bonds("8x79","A","R")
    # bonds("8x79","B","G")
    # bonds("8x79","B","N")
    # bonds("8x79","B","R")

    # 7rg9
    # bonds("7rg9","A","B")

    find.pairs = function(id){
        html = read_html(paste("https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=",id,"&template=interfaces.html&c=999", sep = ""))
        menuClass = html %>%  html_nodes(xpath = sprintf("//*[@class='%s']", "menuClass")) 
        extract_chain = function(x){
        node = html_nodes(x,"img")
        attr = html_attr(node, "src")
    }
    imgs = menuClass %>% sapply(FUN = extract_chain) %>% as.data.frame() %>% mapply(FUN = str_extract, pattern = "chain[a-zA-Z]") %>% as.data.frame() %>% mapply(FUN = gsub, pattern = "chain", replacement = "") %>% t() %>% as.data.frame() %>% rename(chain1 = V1, chain2 = V2) %>% unique()
    rownames(imgs) = NULL
    return(imgs)
    }

find.pairs("6CSX")


    bonds.from.id = function(id){
        id_str = id
        pairs = find.pairs(id)
        bonds = mapply(FUN = bonds, pairs$chain1, pairs$chain2, id = id_str)
    }




    # bonds.from.id("1xhm")
    # bonds.from.id("6n85")
    # bonds.from.id("8kh4")
    # bonds.from.id("8kgk")
    # bonds.from.id("7f8v")
    # bonds.from.id("7yk6")
    # bonds.from.id("7xxi")
    # bonds.from.id("7wvv")
    # bonds.from.id("7wvx")
    # bonds.from.id("6c9h")
    # bonds.from.id("8thk")
    # bonds.from.id("8thl")
    # bonds.from.id("7wvy")
    # bonds.from.id("2ode")
    # bonds.from.id("7x6i")
    # bonds.from.id("8gvx")
    # bonds.from.id("2v4z")
    # bonds.from.id("7e9h")
    # bonds.from.id("8jd6")
    # bonds.from.id("8k9k")
    # bonds.from.id("8szh")
    # bonds.from.id("9avl")
    # bonds.from.id("7t11")
    # bonds.from.id("7t10")
    # bonds.from.id("7ra3")
    # bonds.from.id("8szi")
    # bonds.from.id("7kh0")b
    # bonds.from.id("7ys6")

# gnas_list = c("7E5E", "7BPH", "7YS6", "8X79", "7VUH", "8POK", "8X2K", "8X7A", "7F55", "7F58", "7PIU", "7VUI", "7VUJ", "7XOV", "7Y35", "8GW8", "7XOU", "7PIV", "7RBT", "7RG9", "7RGP", "8EL7", "7F53", "7F54", "7TMW", "7F0T", "7F1O", "7F1Z", "7F23", "7F24", "7RA3", "7X2C", "7Y36", "6XOX", "7KH0")
# this won't work for some reason: "8YUT"
# does not appear to exist in the ebi PDB-- fine 
# gnas_list = c("8ITM", "7X2D", "7X2F", "7X8R", "7X8S", "7XJH", "7XJI")
# 8ZR5 is also not present
# neither is 8ZRK
# gnas_list = c("7FIM", "7WUQ", "7XKD", "7XKF", "7Y3G", "7YP7", "8HMP", "8ITL", "8JLO", "8JLZ", "8PM2", "8TB0", "8U8F", "8W88", "8I2G", "6EG8", "6LI3", "7BW0", "7DH5", "7EVM", "7FIY", "7WCM", "7WCN", "7XTQ", "8HIX", "8HJ0", "8HJ1", "8HJ2", "8HMV", "8IQ4", "8IQ6", "8IZB", "8JHB", "8JIT", "8JIU", "8JPP", "8KGK", "8KH4", "8KH5", "8SMV", "8W8Q", "8WA3", "8WG7", "9BKK", "6ORV", "7AUE", "7CFM", "7CX2")
# gnas_list = c("7CX3", "7CX4", "7DUR", "7F16", "7LCI", "7V35", "7VBH", "7VBI", "7XT9", "7XTB", "7XTC", "8F76")
# gnas_list = c("8FU6", "6B3J", "6WHC", "6WPW", "6X19", "6X1A", "7CFN", "7D3S", "7FIN", "7LLL", "7MBX", "7RTB", "7S1M", "7S3I", "7V9M", "7VQX", "7WBJ", "7XT8", "7XY6", "7XY7", "8FLS", "8FLT", "8FLU", "8GD9", "8GDA", "8GDB", "8JHI", "8JIR", "5UZ7", "5VAI", "6NIY", "6WI9", "6WZG", "7D7M", "7DTY", "7EZK", "7F4I", "7LLY", "7MBY", "7T9N", "7VVJ", "8E3Z", "8FLR", "8GY7", "8HDO", "8HDP", "8IOD")
    gnas_list = c("8W8A", "8ZH8", "6VN7", "7CZ5", "7DHI", "7DHR", "7DUQ", "7F4D", "7F4F", "7F4H", "7RMI", "7TYN", "8E3Y", "8K8J", "8U26", "9ASB", "9AVG", "9BKJ", "3SN6", "6LMK", "6P9Y", "6VCB", "6X18", "7BZ2", "7KI0", "7KI1", "7RMH", "7T9I", "7TYI", "7TYL", "7UTZ", "7VVK", "7VVL", "7VVM", "7VVN", "7VVO", "8E3X", "8F0J", "8GDZ", "8GE1", "8GE2", "8GE3", "8GE4", "8GE5", "8GE6", "8GE7", "8GE8", "8GE9", "8GEA", "8GEB", "8GEC", "8GED", "8GEE", "8GEF", "8GEG", "8GEH", "8GEI", "8GEJ", "8GFV", "8GFW", "8GFX", "8GFY", "8GFZ", "8GG0", "8GG1", "8GG2", "8GG3", "8GG4", "8GG5", "8GG6", "8GG7", "8GG8", "8GG9", "8GGA", "8GGB", "8GGC", "8GGE", "8GGF", "8INR", "8J9N", "8UNL", "8UNM", "8UNN", "8UNO", "8UNP", "8UNQ", "8UNR", "8UNS", "8UNT", "8UNU", "8UNV", "8UNW", "8UNX", "8UNY", "8UNZ", "8UO0", "8UO1", "8UO2", "8UO3", "8UO4", "6GDG", "6NI3", "6UVA", "7BB6", "7BB7", "7CKW", "7CKX", "7CKY", "7CKZ", "7CRH", "7JOZ", "7JV5", "7JVP", "7P00", "7P02", "7TYO", "7TYW", "8FLQ", "8IOC", "9AXF", "6E3Y", "6PB1", "7TYF", "7TYY", "9AUC", "6UUN", "6UUS", "7TYH", "7TYX", "7TZF", "8F0K", "8F2A", "8F2B", "6P9X", "6PB0", "8THK", "8THL", "8UHB") %>% 
        sapply(tolower) %>%
        as.vector()
    # sapply(gnas_list, bonds.from.id)


    # install.packages("pbapply")
    library(pbapply)

    pbsapply(gnas_list, FUN = bonds.from.id)


