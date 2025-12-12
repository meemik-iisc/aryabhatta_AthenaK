# main.py
from menu import select_analysis, build_params_for
from analysis import run

def main():
    analysis = select_analysis()
    print(analysis)
    if analysis=="quit":
        print("Bye!"); return
    params = build_params_for(analysis)
    run(analysis, params)

if __name__=="__main__":
    main()
