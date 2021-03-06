为了方便初学者的同学规范使用pFind 3 软件，董老师计划写一篇使用pFind 3 protocol文章。其基本宗旨是 1）保证新手同学依照文章的步骤可以独立使用软件，遇到问题可以找到对应解决办法。2）涵盖尽量多的使用场景，如开放搜索，限定搜索，定量以及结果解读。
目前，光灿以及整理好了一版开放搜索的step-by-step流程。其实不同使用场景，其主要区别在开放还是限定，以及在限定搜索模式下，如何选择修饰。接下来，我们将要扩充以下内容：
1）限定式搜索
限定式搜索的修饰选择有常见的两种方式，1.根据以往经验进行设置，如C+57作为固定修饰，M+16作为可变修饰；2. 可以根据open search报告的结果作为参考，将丰度较高的修饰作为可变修饰，这里产生了新的问题，选多少个修饰？我们知道可变修饰数目不是越多越好，修饰设置多，一定程度可以提高谱图解析率，但同时也会带来假阳性结果。那么如何给用户提供一个参考的修饰数目？我和光灿打算使用N15标记的数据，根据pfind开放搜索结果分别设置排名靠前的2，3，4，5，6，7，8种修饰来搜索，根据搜索结果报告的肽段数和蛋白数以及谱图层次的非数比例寻找一个平衡的选择。
限定式搜索里一个还有个经常使用的场景是分馏为多个组分的样品数据合在一起搜索。每个raw可能有其特性，不同raw的修饰丰度差异大，我们用怎样的标准去选择修饰使其能够涵盖所有raw的高丰度修饰？
2）定量
pQuant定量的功能非常强大，目前在pFind里也有很好的集成。pFind里面的定量功能包含了：N15，SILAC，labelfree，以及report-intensity（TMT)。这些都是常用的定量方法，有必要在文章中一一展示。然而后面3种我们几乎没有做过（需要数据支持），对这些数据的理解也不深，所以要写也只是把整个过程走一遍。

