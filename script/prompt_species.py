from openai import OpenAI

with open("api_key.txt", "r", encoding="utf-8") as f:
    api_key = f.read().strip()
client = OpenAI(api_key=api_key)

species_name = input("Enter species name: ")

with open("prompt.txt", "r", encoding="utf-8") as file:
    prompt_template = file.read()

prompt = prompt_template.replace("{species}", species_name)

response = client.chat.completions.create(
    model="gpt-4o-mini",
    messages= [{"role": "user", "content": prompt}],
)

print("\nModel Response:\n")
print(response.choices[0].message.content)
