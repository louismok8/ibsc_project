from models.template import BiopsyTemplate

def main():
    template = BiopsyTemplate(
        grid_size=19,
        spacing=5.0,
        origin=(0.0, 0.0, 0.0),
        direction=(0.0, 0.0, 1.0),
    )

    print(f"Total grid holes: {len(template)}")
    print("First hole:", template.holes[0])
    print("Center hole:", template.holes[len(template)//2])

if __name__ == "__main__":
    main()
