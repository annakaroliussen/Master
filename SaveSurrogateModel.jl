using Serialization

# File where all surrogate models are stored
const MODEL_FILE = "/home/annkar/NewCode/Data/Dictionary/surrogate_models_oneInput.jls"
const MODEL_FILE1 = "/home/annkar/NewCode/Data/Dictionary/surrogate_models_twoInput.jls"

# =========================================================
# Save surrogate model
# =========================================================

function SaveModel_oneInput(key, beta, c, SSE)

    # Load existing models
    if isfile(MODEL_FILE)
        models = deserialize(MODEL_FILE)
    else
        models = Dict()
    end

    # Save new model
    models[key] = Dict(
        "beta" => beta,
        "c" => c,
        "SSE" => SSE
    )

    # Save dictionary to file
    serialize(MODEL_FILE, models)

    println("Model saved: ", key)

end


function SaveModel_TwoInput(key, beta, c, SSE)

    # Load existing models
    if isfile(MODEL_FILE1)
        models = deserialize(MODEL_FILE1)
    else
        models = Dict()
    end

    # Save new model
    models[key] = Dict(
        "beta" => beta,
        "c" => c,
        "SSE" => SSE
    )

    # Save dictionary to file
    serialize(MODEL_FILE1, models)

    println("Model saved: ", key)

end


# =========================================================
# Load surrogate model
# =========================================================

function LoadModel_OneInput(key)

    models = deserialize(MODEL_FILE)

    return models[key]

end
function LoadModel_TwoInput(key)

    models = deserialize(MODEL_FILE1)

    return models[key]

end


# =========================================================
# Show all saved models
# =========================================================

function ShowModels_OneInput()

    if !isfile(MODEL_FILE)
        println("No models saved yet.")
        return
    end

    models = deserialize(MODEL_FILE)

    println("Saved surrogate models:")

    for key in keys(models)
        println(" - ", key)
    end

end

function ShowModels_TwoInput()

    if !isfile(MODEL_FILE1)
        println("No models saved yet.")
        return
    end

    models = deserialize(MODEL_FILE1)

    println("Saved surrogate models:")

    for key in keys(models)
        println(" - ", key)
    end

end