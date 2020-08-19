name: CompatHelper
on:
  schedule:
    - cron: '00 00 * * *'
  workflow_dispatch:
jobs:
  CompatHelper:
    runs-on: ubuntu-latest
    steps:
      - name: Pkg.add("CompatHelper")
        run: julia -e 'using Pkg; Pkg.add("CompatHelper")'
      - name: CompatHelper.main()
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          COMPATHELPER_PRIV: ${{ secrets.COMPATHELPER_PRIV }}  # optional
        run: julia -e 'using CompatHelper; CompatHelper.main()'
