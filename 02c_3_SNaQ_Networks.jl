# Script for Phylogenetic Network Analysis Using PhyloNetworks
# -----------------------------------------------------------

# 1. Install necessary packages (uncomment if needed)
# Pkg.add(PackageSpec(name="PhyloNetworks", version="0.16.4"))
# Pkg.add("PhyloPlots")
# Pkg.add("RCall")      
# Pkg.add("CSV")        
# Pkg.add("DataFrames") 
# Pkg.add("StatsModels")

# 2. Load required libraries
using PhyloNetworks
using PhyloPlots
using CSV
using DataFrames
using RCall
using QuartetNetworkGoodnessFit

# 3. Set parameters
group = " "  # outfile prefix, also the name of the directory with all the input files
out = " "  # outgroup

# 4. Setup directory structure
BASE_DIR = joinpath(pwd(), "PhyloReconcile")
MSC_RESULTS = joinpath(BASE_DIR, "09_Network_Analyses")
group_dir = joinpath(MSC_RESULTS, group)

# 5. Define input and output file paths
tree_file = joinpath(group_dir, string(group, ".tre"))
cf_file = joinpath(group_dir, string(group, "_individualsCFs.csv"))
INPUTS = joinpath(BASE_DIR, "01_initial_data", "input_files", "SNaQ_NANUQ")
mapping_file = joinpath(INPUTS, string(group, "_map.csv"))
sptable_file = joinpath(group_dir, string(group, "CFtable_species.csv"))
sptable_norep_file = joinpath(group_dir, string(group, "CFtable_species_norep.csv"))
results_file = joinpath(group_dir, string("results", group, ".csv"))
astraltree = readTopology(tree_file)

# 7. Make a species-level CF file, write to folder
df_sp = mapAllelesCFtable(mapping_file, cf_file)
d_sp = readTableCF!(df_sp)
df_sp = writeTableCF(d_sp)
CSV.write(sptable_file, df_sp)

# 8. Run SNaQ with increasing hybridization levels
# h=0 (no hybridization)
numretic = 0
output_file_0 = joinpath(group_dir, string("net", numretic, group))
net0 = snaq!(astraltree, d_sp, hmax=0, filename=output_file_0, seed=2345)

# h=1
numretic = 1
output_file_1 = joinpath(group_dir, string("net", numretic, group))
net1 = snaq!(net0, d_sp, hmax=1, filename=output_file_1, seed=3456)

# h=2
numretic = 2
output_file_2 = joinpath(group_dir, string("net", numretic, group))
net2 = snaq!(net1, d_sp, hmax=2, filename=output_file_2, seed=4567)

# h=3
numretic = 3
output_file_3 = joinpath(group_dir, string("net", numretic, group))
net3 = snaq!(net2, d_sp, hmax=3, filename=output_file_3, seed=5678)

# h=4
numretic = 4
output_file_4 = joinpath(group_dir, string("net", numretic, group))
net4 = snaq!(net3, d_sp, hmax=4, filename=output_file_4, seed=6789)

# h=5
numretic = 5
output_file_5 = joinpath(group_dir, string("net", numretic, group))
net5 = snaq!(net4, d_sp, hmax=5, filename=output_file_5, seed=7891)

# 9. Process species-level file and remove repeats (1 representative per species)
df_sp = DataFrame(CSV.File(sptable_file))

function hasrep(row)
    occursin(r"__2$", row[:t1]) || occursin(r"__2$", row[:t2]) ||
    occursin(r"__2$", row[:t3]) || occursin(r"__2$", row[:t4])
end

df_sp_reduced = filter(!hasrep, df_sp)
CSV.write(sptable_norep_file, df_sp_reduced)

d_sp_reduced = readTableCF(df_sp_reduced)  # From PhyloNetworks
d_sp = d_sp_reduced

qCF = DataFrame(CSV.File(sptable_norep_file))

# 10. Read in tree and calculate maxq pseudolik and goodness of fit
net0_file = joinpath(group_dir, string("net0", group, ".out"))
net0 = readTopology(net0_file)
net0alt = topologyMaxQPseudolik!(net0, d_sp)
ressnaq = quarnetGoFtest!(net0, qCF, false; seed=234, nsim=200)
ressnaq1 = quarnetGoFtest!(net0alt, qCF, false; seed=234, nsim=200)
x = writeTopology(net0alt)

# 11. Create dataframe for the results
dffull = DataFrame(
    network = "0",
    numhybrids = 0,
    snaq = net0.loglik, 
    maxq = net0alt.loglik, 
    gof = ressnaq[1], 
    gofalt = ressnaq1[1], 
    topology = "x"
)

# 12. Process all possible reticulation searches (1 to 5)
for j in 1:5
    network_file = joinpath(group_dir, string("net", j, group, ".networks"))
    
    if isfile(network_file)  # if there is a reticulation search
        # Read in the topologies
        networks = readMultiTopology(network_file)
        
        # Get the loglikelihoods - FIXED SED COMMAND
        # The previous command was causing problems with the '\1' capture group
        try
            # Alternative approach to extract loglik values
            file_content = read(network_file, String)
            # Find all occurrences of loglik values
            loglik_pattern = r"with -loglik ([0-9]+\.[0-9]+)"
            loglik_matches = collect(eachmatch(loglik_pattern, file_content))
            scores1 = [parse(Float64, m.captures[1]) for m in loglik_matches]
            
            numnets = length(networks)
            @rput numnets
            
            # Create individual SVG files for each network
            for i in eachindex(networks)
                # Get rid of hybrids <10% contribution, calculate the new maximized pseudologliklihood, and new goodness of fit value
                print(string(j, ".", i, "\n"))
                deleteHybridThreshold!(networks[i], 0.1)
                x = writeTopology(networks[i])
                
                net1alt = try
                    topologyMaxQPseudolik!(networks[i], d_sp)
                catch
                    0
                end
                
                loglik = try
                    net1alt.loglik
                catch
                    0
                end
                
                ressnaq = try
                    quarnetGoFtest!(networks[i], qCF, false; seed=234, nsim=200)
                catch
                    0
                end
                
                resalt = try
                    quarnetGoFtest!(net1alt, qCF, false; seed=234, nsim=200)
                catch
                    0
                end
                
                ressnaq1 = try
                    ressnaq[1]
                catch
                    0
                end
                
                resalt1 = try
                    resalt[1]
                catch
                    0
                end
                
                # Put all of that info in the data frame
                push!(dffull, (string(j, ".", i), networks[i].numHybrids, scores1[i], loglik, ressnaq1, resalt1, x))
                
                # For plotting - create individual SVG for each network
                @rput group
                @rput j
                @rput i
                
                svg_file = joinpath(group_dir, string(group, j, "_", i, "maxq.svg"))
                @rput svg_file
                R"svg(svg_file, width=7, height=4)"
                R"par(mar=c(0,3,3,0))"  # Set margins
                R"par(cex=0.75)"        # Set font size
                
                # For plotting
                try
                    rootatnode!(networks[i], out)
                catch e
                    0
                end
                
                try
                    rootatnode!(net1alt, out)
                catch e
                    0
                end
                
                try
                    plot(net1alt, :R, showGamma=true)
                    num = string(j, ".", i, " loglik=", scores1[i], " maxqLL=", loglik, " p=", ressnaq1, " maxqp=", resalt1)
                    @rput num
                    R"mtext"(num, cex=0.6)
                catch e
                    plot(networks[i], :R, showGamma=true)
                    plot(net0, :R)
                    num = string(j, ".", i, " loglik=", scores1[i], " maxqLL=", loglik, " p=", ressnaq1, " maxqp=", resalt1)
                    @rput num
                    R"mtext"(num, cex=0.6)
                end
                
                # Close the SVG file after each plot
                R"dev.off()"
            end
        catch e
            println("Error processing network file: ", network_file)
            println("Error: ", e)
        end
    end
    
    # Write to file as it updates
    CSV.write(results_file, dffull)
end

println("Analysis completed for group: ", group)
println("Results saved to: ", results_file)
